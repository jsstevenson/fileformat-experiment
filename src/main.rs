use std::path::PathBuf;
pub mod vcf;
pub mod parse;
use tokio;

#[tokio::main]
async fn main() {
    println!("hello world");

    let path = PathBuf::from(
        r"/Users/jss009/code/vrs_anvil_toolkit/u08_release_data/gregor_consortium_u06_sorted_chrY_V2_VT_VEP_VRS.vcf.gz",
    );
    let _ = vcf::stream_from_vcf(path).await;
    ()
}
