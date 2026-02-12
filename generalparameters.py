# here we have the general parameters for a run

INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/LONGAPPROACH"
N_FACTORS_LIST = list(range(10, 31, 2))

VIEWS_CONFIG = {
    "Metabolomics": {"path": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Lipidomics": {"path": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Microbiota": {"path": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv", "likelihood": "Normal", "scale": True},
    "SV_DEL": {"path": f"{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_DUP": {"path": f"{INPUT_DIR}/SV_DUP_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_INS": {"path": f"{INPUT_DIR}/SV_INS_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_INV": {"path": f"{INPUT_DIR}/SV_INV_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_TRA": {"path": f"{INPUT_DIR}/SV_TRA_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "VC_11": {"path": f"{INPUT_DIR}/MPA_GT_1_1.tsv", "likelihood": "Bernoulli", "scale": False},
    "VC_12": {"path": f"{INPUT_DIR}/MPA_GT_1_2.tsv", "likelihood": "Bernoulli", "scale": False},
    "EXPRESSION": {"path": f"{INPUT_DIR}/tpmexpression.tsv", "likelihood": "Normal", "scale": True}
}
