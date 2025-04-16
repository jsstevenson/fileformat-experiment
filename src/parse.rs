use noodles_vcf::{
    self as vcf,
    variant::record::info::{self, field::Value as InfoValue},
};

#[derive(Debug)]
pub struct AlleleRow {
    file_uri: String,
    chr: String,
    pos: u32, // good enough for humans
    vrs_id: String,
    vrs_start: u32,
    vrs_end: u32,
    vrs_state: String
}

#[derive(Debug)]
pub enum VcfParseError {
    Error(String)
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

/// Extracts an array of strings from a VCF record's info field.
///
/// This function attempts to retrieve the info field specified by `field` from the given VCF record info
/// using the provided `header`. It expects the field to be an array of strings, and if successful,
/// returns a vector of `String` values. If the field is not found, not an array, or is not stored as strings,
/// an error is returned.
///
/// # Arguments
///
/// * `info` - The VCF record info object.
/// * `header` - A reference to the VCF header used to interpret the info field.
/// * `field` - The field identifier as a `VrsVcfField`.
///
/// # Returns
///
/// A `Result` containing a vector of strings if the info field is successfully extracted and converted;
/// otherwise, returns a `PyErr`.
///
/// # Errors
///
/// Returns an error if:
/// - The specified field is not present in the info.
/// - The field is present but not stored as an array of strings.
/// - The array variant does not match the expected type.
pub fn get_str_info_field(
    info: vcf::record::Info,
    header: &vcf::Header,
    field: VrsVcfFieldName,
) -> Result<Vec<String>, VcfParseError> {
    if let Some(Ok(Some(InfoValue::Array(ids_array)))) = info.get(header, field.as_str()) {
        if let info::field::value::Array::String(array_elements) = ids_array {
            let vec = array_elements
                .iter()
                .map(|cow_str| cow_str.unwrap().unwrap_or_default().to_string())
                .collect();
            return Ok(vec);
        } else {
            Err(VcfParseError::Error("expected string array variant".to_string()))
        }
    } else {
        Err(VcfParseError::Error("Expected Array variant".to_string()))
    }
}

/// Extracts an integer array from a VCF record's info field.
///
/// This function attempts to retrieve the info field specified by `field` from the VCF record `info`
/// using the provided `header`. It expects the field to be an array stored either as integers or as strings.
///
/// - If the array is stored as integers, the function unwraps each element and returns a `Vec<i32>`.
/// - If the array is stored as strings, each string is parsed into an `i32`.
///
/// The VRS-Python annotator now processes numeric annotation fields as Integers, but it didn't
/// always, so the fallback is there to manage legacy VCF outputs.
///
/// # Arguments
///
/// * `info` - The VCF record info object.
/// * `header` - A reference to the VCF header, which is used to interpret the info field.
/// * `field` - The field identifier as a `VrsVcfField`.
///
/// # Returns
///
/// Returns `Ok(Vec<i32>)` if the info field is successfully extracted and parsed. Otherwise, an error is returned:
/// - If the field is not present or not an array.
/// - If the array is not in an expected format (i.e. neither integer nor string array).
/// - If parsing a string into an integer fails.
///
/// # Errors
///
/// This function returns a `PyErr` if the expected array variant is not found or if any parsing error occurs.
pub fn get_int_info_field(
    info: vcf::record::Info,
    header: &vcf::Header,
    field: VrsVcfFieldName,
) -> Result<Vec<i32>, VcfParseError> {
    if let Some(Ok(Some(InfoValue::Array(ids_array)))) = info.get(header, field.as_str()) {
        if let info::field::value::Array::Integer(array_elements) = ids_array {
            let vec = array_elements
                .iter()
                .map(|i| i.unwrap().unwrap_or_default())
                .collect();
            return Ok(vec);
        } else {
            if let info::field::value::Array::String(array_elements) = ids_array {
                let vec = array_elements
                    .iter()
                    .map(|cow_str| cow_str.unwrap().unwrap_or_default().to_string())
                    .map(|s| {
                        s.parse::<i32>()
                            .map_err(|_e| VcfParseError::Error("Unable to parse int".to_string()))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                return Ok(vec);
            }
            Err(VcfParseError::Error(
                "Expected array variable of strings or integers".to_string(),
            ))
        }
    } else {
        Err(VcfParseError::Error("Expected Array variant".to_string()))
    }
}
