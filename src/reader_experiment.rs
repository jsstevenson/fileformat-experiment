use futures::TryStreamExt;
use noodles_bgzf::r#async::Reader as BgzfReader;
use noodles_vcf::{
    self as vcf,
    r#async::io::Reader as VcfReader,
    variant::record::info::{self, field::Value as InfoValue},
};
use std::path::PathBuf;
use std::time::Instant;
use tokio::{
    fs::File as TkFile,
    io::{AsyncBufRead, BufReader},
};

enum Error {
    MyErr,
}

fn noodles_get_vrs_ids(
    info: vcf::record::Info,
    header: &vcf::Header,
) -> Result<Vec<String>, Error> {
    if let Some(Ok(Some(InfoValue::Array(array)))) = info.get(header, "VRS_Allele_IDs") {
        if let info::field::value::Array::String(array_elements) = array {
            let vec = array_elements
                .iter()
                .map(|cow_str| cow_str.unwrap().unwrap_or_default().to_string())
                .collect();
            return Ok(vec);
        } else {
            Err(Error::MyErr)
        }
    } else {
        Err(Error::MyErr)
    }
}

async fn noodles_get_reader(
    vcf_path: PathBuf,
) -> Result<VcfReader<Box<dyn tokio::io::AsyncBufRead + Unpin + Send>>, Error> {
    let file = TkFile::open(vcf_path.clone())
        .await
        .map_err(|_| Error::MyErr)?;
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
        _ => Err(Error::MyErr),
    }
}

#[derive(Debug)]
pub struct DbRow {
    pub vrs_id: String,
    pub chr: String,
    pub pos: i64,
}

async fn noodles_load_vcf(vcf_path: PathBuf) -> Result<(), Error> {
    let start = Instant::now();

    let mut reader = noodles_get_reader(vcf_path).await?;
    let header = reader.read_header().await.map_err(|_| Error::MyErr)?;

    let mut records = reader.records();

    while let Some(record) = records.try_next().await.map_err(|_| Error::MyErr)? {
        let vrs_ids = noodles_get_vrs_ids(record.info(), &header)?;
        let chrom = record.reference_sequence_name();
        let pos = record.variant_start().unwrap().unwrap().get();

        for vrs_id in vrs_ids {
            let _ = DbRow {
                vrs_id: vrs_id
                    .strip_prefix("ga4gh:VA.")
                    .unwrap_or(&vrs_id)
                    .to_string(),
                chr: chrom.to_string(),
                pos: pos.try_into().unwrap(),
            };
        }
    }

    let duration = start.elapsed();
    println!("Time taken: {:?}", duration);
    Ok(())
}

use std::fs::File;
use std::io::{self, BufRead, BufReader as StdBufReader};

fn naive_load_vcf(vcf_path: PathBuf) -> Result<(), Error> {
    let start = Instant::now();
    let file = File::open(vcf_path).unwrap();
    let reader = StdBufReader::new(file);
    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') {
            continue; // Skip header lines (## and #CHROM)
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue; // Skip lines that don't have at least 8 columns
        }

        let chr = fields[0].to_string();
        let pos = match fields[1].parse::<i64>() {
            Ok(n) => n,
            Err(_) => continue, // Skip if position isn't an integer
        };

        let mut kv_pairs = fields[7].split(';');
        let vrs_ids_value = kv_pairs.find_map(|kv| {
            let mut split = kv.splitn(2, '=');
            match (split.next(), split.next()) {
                (Some(key), Some(value)) if key == "VRS_Allele_IDs" => Some(value),
                _ => None,
            }
        });

        if let Some(vrs_ids_str) = vrs_ids_value {
            for vrs_id in vrs_ids_str.split(',') {
                let _ = DbRow {
                    vrs_id: vrs_id.trim().to_string(),
                    chr: chr.clone(),
                    pos,
                };
            }
        }
    }

    let duration = start.elapsed();
    println!("Time taken: {:?}", duration);
    Ok(())
}

#[tokio::main()]
async fn main() {
    let path = PathBuf::from(
        r"/Users/jss009/code/vrs_anvil_toolkit/u08_release_data/gregor_consortium_u06_sorted_chrY_V2_VT_VEP_VRS.vcf.gz",
    );
    let _ = noodles_load_vcf(path.clone()).await;
    let _ = naive_load_vcf(path.clone()); //  TODO how to gunzip
}
