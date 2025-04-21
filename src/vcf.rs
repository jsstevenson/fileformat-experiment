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
    ParseFailure(String),
    NullField,
    TmpErr, // placeholder, basically
}

/// Represents the different kind of supported VRS variations
///
/// Currently, the VCF annotator can only translate to alleles, so other variation
/// types supported by VCF get dropped.
enum VariationType {
    Allele,
    // TODO others
}

impl VariationType {
    /// Returns an ID for each variation type.
    ///
    /// When we output to the compressed file, this lets us e.g. treat "VA." as "1".
    ///
    /// # Errors
    ///
    /// Returns an error if the mapping fails (This is ~impossible)
    fn to_id(&self) -> Result<u8, ()> {
        match self {
            VariationType::Allele => Ok(1),
        }
    }
}

/// Contains a single set of VRS attributes grabbed from an INFO field.
#[derive(Debug)]
struct VrsAlleleAttrs {
    vrs_id: String,
    vrs_start: i32,
    vrs_end: i32,
    vrs_state: String,
}

impl VrsAlleleAttrs {
    /// Convert VRS ID to compressed form
    ///
    /// strip ga4gh: if it's there, change the type value to a shortened ID)
    ///
    /// # Errors
    ///
    /// If unrecognized prefix is encountered (this should be impossible)
    fn vrs_id_to_vrsix(&self) -> Result<String, VcfError> {
        let no_namespace = self.vrs_id.strip_prefix("ga4gh:").unwrap_or(&self.vrs_id);
        match no_namespace {
            s if s.starts_with("VA.") => {
                let rest = s.strip_prefix("VA.").unwrap();
                let rep = VariationType::Allele.to_id().unwrap();
                Ok(format!("{}{}", rep, rest))
            }
            _ => Err(VcfError::TmpErr),
        }
    }
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
                Ok(None) => Ok("".to_string()),
                Err(_) => Err(VcfError::ParseFailure(
                    "Individual array element failed to parse".to_string(),
                )),
            });
            let collected: Result<Vec<_>, _> = iter.collect();
            collected.map(|vec| vec.into_iter())
        } else {
            Err(VcfError::ParseFailure(
                "Failed to parse as array of strings".to_string(),
            ))
        }
    } else {
        Err(VcfError::ParseFailure(
            "Failed to parse as array".to_string(),
        ))
    }
}

fn get_vrs_pos(
    info: vcf::record::Info,
    header: &vcf::Header,
    field: VrsVcfFieldName,
) -> Result<impl Iterator<Item = i32>, VcfError> {
    if let Some(Ok(Some(InfoValue::Array(array)))) = info.get(header, field.as_str()) {
        match array {
            info::field::value::Array::Integer(array_elements) => {
                let iter = array_elements.iter().map(|res_opt| match res_opt {
                    Ok(Some(num)) => Ok(num),
                    Ok(None) => Err(VcfError::TmpErr), // TODO handle this case
                    Err(_) => Err(VcfError::ParseFailure(
                        "Individual array element failed to parse".to_string(),
                    )),
                });
                let collected: Result<Vec<_>, _> = iter.collect();
                collected.map(|v| v.into_iter())
            }
            // handle old cases where the position columns were strings
            info::field::value::Array::String(array_elements) => {
                let iter = array_elements.iter().map(|res_opt| match res_opt {
                    Ok(Some(cow)) => Ok(cow.to_string().parse::<i32>().unwrap()),
                    Ok(None) => Err(VcfError::TmpErr), // TODO handle this case
                    Err(_) => Err(VcfError::ParseFailure(
                        "Individual array element failed to parse".to_string(),
                    )),
                });
                let collected: Result<Vec<_>, _> = iter.collect();
                collected.map(|vec| vec.into_iter())
            }
            _ => Err(VcfError::ParseFailure(
                "Failed to parse as array of ints".to_string(),
            )),
        }
    } else {
        Err(VcfError::TmpErr)
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
                vrs_start,
                vrs_end,
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

#[derive(Debug)]
struct FileData {
    chrom: String,
    pos: u32,
    uri_id: u8,
    vrs_hash: String, // 1 byte type ID + 32 bytes, ASCII
    vrs_start: i32,
    vrs_end: i32,
    vrs_state: String, // varchar but should be ASCII
}

// File layout
// <header>
// <records -- uri_id, chrom, pos
// <vartype + vrs_id, seek offset>
// <vrs start, seek offset>
// <vrs end, seek offset>
use tokio::fs::{File, OpenOptions};
use tokio::io::AsyncWriteExt;

enum OutfileError {
    General
}


async fn write_data_to_file(out_file: &mut tokio::fs::File, line: String) -> Result<(), OutfileError> {
    out_file.write_all(line.as_bytes())
        .await
        .map_err(|_| OutfileError::General)?;
    out_file.flush()
        .await
        .map_err(|_| OutfileError::General)?;
    Ok(())
}

pub async fn load_vcf(vcf_path: PathBuf, file_uri: Option<String>, output_file: PathBuf) -> Result<(), VcfError> {
    let mut reader = get_reader(vcf_path)
        .await
        .map_err(|_| VcfError::TmpErr)
        .unwrap();
    let header = reader.read_header().await.unwrap();

    let mut records = reader.records();
    let mut out_file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(output_file)
        .await
        .map_err(|_| VcfError::TmpErr)?;
    let mut count = 0;

    let uri_id: u8 = 1;  // TODO figure out how to calculate this

    while let Some(record) = records.try_next().await.map_err(|_| VcfError::TmpErr)? {
        let chrom = record.reference_sequence_name();
        let pos = record.variant_start().unwrap().unwrap().get() as u32;

        let mut stream = record.iter_vrs_attrs(&header).await;
        while let Some(attrs_result) = stream.next().await {
            match attrs_result {
                Ok(attrs) => {
                    let data = FileData {
                        chrom: chrom.to_string(),
                        pos,
                        uri_id,
                        vrs_hash: attrs.vrs_id_to_vrsix().unwrap(),
                        vrs_start: attrs.vrs_start,
                        vrs_end: attrs.vrs_end,
                        vrs_state: attrs.vrs_state
                    };
                    let line = format!("{}-{}-{}\n", data.chrom, data.pos, data.uri_id);
                    let _ = write_data_to_file(&mut out_file, line).await;
                    let line = format!("{}{}\n", data.vrs_hash, count);
                    let _ = write_data_to_file(&mut out_file, line).await;
                    let line = format!("{}-{}\n", data.vrs_start, count);
                    let _ = write_data_to_file(&mut out_file, line).await;
                    let line = format!("{}-{}\n", data.vrs_end, count);
                    let _ = write_data_to_file(&mut out_file, line).await;
                    count += 1;
                }
                Err(attrs) => eprintln!("{:?}", attrs),
            }
        }
    }
    Ok(())
}
