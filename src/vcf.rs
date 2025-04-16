use async_trait::async_trait;
use futures::{stream, stream::StreamExt, TryStreamExt};
use itertools::multizip;
use noodles_bgzf::r#async::Reader as BgzfReader;
use noodles_vcf::{
    self as vcf,
    r#async::io::Reader as VcfReader,
    variant::record::info::{self, field::Value as InfoValue},
    Record,
};
use std::path::PathBuf;
use tokio::{
    fs::File as TkFile,
    io::{AsyncBufRead, BufReader},
};

#[derive(Debug)]
pub enum VcfError {
    UnsupportedFiletype,
    MyErr, // placeholder, basically
}

// attrs to stream from a VCF row
#[derive(Debug)]
struct VrsAlleleAttrs {
    vrs_id: String,
    vrs_start: u64,
    vrs_end: u64,
    vrs_state: String,
}

#[derive(Debug)]
pub enum VrsVcfFieldName {
    VrsAlleleIds,
    VrsStarts,
    VrsEnds,
    VrsStates,
}

impl VrsVcfFieldName {
    fn as_str(&self) -> &'static str {
        match self {
            VrsVcfFieldName::VrsAlleleIds => "VRS_Allele_IDs",
            VrsVcfFieldName::VrsStarts => "VRS_Starts",
            VrsVcfFieldName::VrsEnds => "VRS_Ends",
            VrsVcfFieldName::VrsStates => "VRS_States",
        }
    }
}

fn get_vrs_str_field(
    info: vcf::record::Info,
    header: &vcf::Header,
    field: VrsVcfFieldName,
) -> Result<impl Iterator<Item = String>, VcfError> {
    if let Some(Ok(Some(InfoValue::Array(array)))) = info.get(header, field.as_str()) {
        if let info::field::value::Array::String(array_elements) = array {
            let iter = array_elements.iter().map(|res_opt| match res_opt {
                Ok(Some(cow)) => Ok(cow.to_string()),
                Ok(None) => Err(VcfError::MyErr),
                Err(_) => Err(VcfError::MyErr),
            });
            let collected: Result<Vec<_>, _> = iter.collect();
            collected.map(|vec| vec.into_iter())
        } else {
            Err(VcfError::MyErr)
        }
    } else {
        Err(VcfError::MyErr)
    }
}

fn get_vrs_pos(
    info: vcf::record::Info,
    header: &vcf::Header,
    field: VrsVcfFieldName,
) -> Result<impl Iterator<Item = i32>, VcfError> {
    if let Some(Ok(Some(InfoValue::Array(array)))) = info.get(header, field.as_str()) {
        if let info::field::value::Array::Integer(array_elements) = array {
            let iter = array_elements.iter().map(|res_opt| match res_opt {
                Ok(Some(num)) => Ok(num),
                Ok(None) => Err(VcfError::MyErr),
                Err(_) => Err(VcfError::MyErr),
            });
            let collected: Result<Vec<_>, _> = iter.collect();
            collected.map(|v| v.into_iter())
        } else {
            // TODO handle case where strings need to be coerced
            Err(VcfError::MyErr)
        }
    } else {
        Err(VcfError::MyErr)
    }
}

#[async_trait]
trait InfoFieldTranspose {
    async fn iter_vrs_attrs(
        &self,
        header: &vcf::Header,
    ) -> stream::BoxStream<'_, Result<VrsAlleleAttrs, VcfError>>;
}

#[async_trait]
impl InfoFieldTranspose for Record {
    async fn iter_vrs_attrs(
        &self,
        header: &vcf::Header,
    ) -> stream::BoxStream<'_, Result<VrsAlleleAttrs, VcfError>> {
        let vrs_id_iter =
            match get_vrs_str_field(self.info(), header, VrsVcfFieldName::VrsAlleleIds) {
                Ok(iter) => iter,
                Err(e) => return stream::once(async { Err(e) }).boxed(),
            };
        let vrs_start_iter = match get_vrs_pos(self.info(), header, VrsVcfFieldName::VrsStarts) {
            Ok(iter) => iter,
            Err(e) => return stream::once(async { Err(e) }).boxed(),
        };
        let vrs_end_iter = match get_vrs_pos(self.info(), header, VrsVcfFieldName::VrsEnds) {
            Ok(iter) => iter,
            Err(e) => return stream::once(async { Err(e) }).boxed(),
        };
        let vrs_state_iter =
            match get_vrs_str_field(self.info(), header, VrsVcfFieldName::VrsStates) {
                Ok(iter) => iter,
                Err(e) => return stream::once(async { Err(e) }).boxed(),
            };
        let stream = stream::iter(multizip((
            vrs_id_iter,
            vrs_start_iter,
            vrs_end_iter,
            vrs_state_iter,
        )))
        .map(|(vrs_id, vrs_start, vrs_end, vrs_state)| {
            Ok(VrsAlleleAttrs {
                vrs_id,
                vrs_start: vrs_start as u64,
                vrs_end: vrs_end as u64,
                vrs_state,
            })
        });
        stream.boxed()
    }
}

async fn get_reader(
    vcf_path: PathBuf,
) -> Result<VcfReader<Box<dyn tokio::io::AsyncBufRead + Unpin + Send>>, VcfError> {
    let file = TkFile::open(vcf_path.clone()).await.unwrap();
    let ext = vcf_path.extension().and_then(|ext| ext.to_str());
    match ext {
        Some("gz") => {
            let reader = Box::new(BgzfReader::new(file)) as Box<dyn AsyncBufRead + Unpin + Send>;
            Ok(VcfReader::new(reader))
        }
        Some("vcf") => {
            let reader = Box::new(BufReader::new(file)) as Box<dyn AsyncBufRead + Unpin + Send>;
            Ok(VcfReader::new(reader))
        }
        _ => Err(VcfError::UnsupportedFiletype),
    }
}

pub async fn load_vcf(vcf_path: PathBuf) -> Result<(), VcfError> {
    let mut reader = get_reader(vcf_path)
        .await
        .map_err(|_| VcfError::MyErr)
        .unwrap();
    let header = reader.read_header().await.unwrap();

    let mut records = reader.records();

    while let Some(record) = records.try_next().await.map_err(|_| VcfError::MyErr)? {
        let chrom = record.reference_sequence_name();
        let pos = record.variant_start().unwrap().unwrap().get();

        let mut stream = record.iter_vrs_attrs(&header).await;
        while let Some(attrs) = stream.next().await {
            match attrs {
                Ok(attrs) => (),
                Err(attrs) => eprintln!("{:?}", attrs)
            }
        }


    }

    Ok(())
}
