import pandas as pd
import numpy as npls
import warnings
import gc
import os
try:
    import plotnine as pn
    HAVE_PLOTNINE = True
except Exception:
    pn = None
    HAVE_PLOTNINE = False
    print("Warning: plotnine not available; faceted factor plot will use matplotlib fallback.")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    import seaborn as sns
    HAVE_SEABORN = True
except Exception:
    sns = None
    HAVE_SEABORN = False
    print("Warning: seaborn not available; using matplotlib-only plotting fallback.")

# MOFAFLEX imports
import mofaflex as mfl
import mudata as md
import anndata as ad

# Configure matplotlib for headless environment

warnings.filterwarnings("ignore", category=FutureWarning)


# Helper to save plots returned by different plotting backends (matplotlib or plotnine/ggplot)
def save_plot(plot_obj, filename, dpi=300):
    """Save a plot object to file. Supports matplotlib Figure/Axes and plotnine ggplot objects.

    If plot_obj is a matplotlib Figure or Axes, use savefig.
    If it's a plotnine ggplot, use its save method.
    If it's a list/tuple of figures, try to save the first or iterate.
    """
    try:
        # matplotlib Figure
        if hasattr(plot_obj, 'savefig'):
            plot_obj.savefig(filename, dpi=dpi)
            return True
        # matplotlib Axes -> get the figure
        if hasattr(plot_obj, 'figure') and hasattr(plot_obj.figure, 'savefig'):
            plot_obj.figure.savefig(filename, dpi=dpi)
            return True
        # plotnine / ggplot
        if hasattr(plot_obj, 'save'):
            # plotnine's save takes filename and dpi
            try:
                plot_obj.save(filename, dpi=dpi)
                return True
            except TypeError:
                # some versions accept only filename
                plot_obj.save(filename)
                return True
        # list/tuple of figures
        if isinstance(plot_obj, (list, tuple)) and len(plot_obj) > 0:
            for i, p in enumerate(plot_obj):
                try:
                    save_plot(p, filename if i == 0 else f"{os.path.splitext(filename)[0]}_{i}{os.path.splitext(filename)[1]}", dpi=dpi)
                except Exception:
                    continue
            return True
    except Exception as e:
        print(f"Could not save plot to {filename}: {e}")
    print(f"Warning: unknown plot object type; could not save {filename}")
    return False


def safe_model_save(model, path, **kwargs):
    """Save a MOFAFLEX model, handling API differences between versions.

    Tries model.save(path) first, falls back to model._save(path, **kwargs).
    """
    try:
        if hasattr(model, 'save') and callable(getattr(model, 'save')):
            return model.save(path)
    except Exception:
        pass
    # fallback to _save
    if hasattr(model, '_save') and callable(getattr(model, '_save')):
        return model._save(path, **kwargs)
    raise AttributeError('MOFAFLEX model has no save or _save method')

# Memory optimization functions
def check_memory_usage():
    """Check current memory usage"""
    try:
        import psutil
        memory = psutil.virtual_memory()
        return f"Memory usage: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB used / {memory.total/1e9:.1f}GB total)"
    except ImportError:
        return "Memory monitoring unavailable"

# Provide analysis summary
print(f"\n📈 ANALYSIS OUTPUTS:")
print(f"  - Optimized datasets saved as CSV files")
print(f"  - MOFAFLEX model trained and saved as H5 file")
print(f"  - Configuration saved in mofa_config.json")

print(f"\n🔧 NEXT STEPS:")
print(f"  1. Load the saved model to explore factors and loadings")
print(f"  2. Consider increasing factors (to 10-12) for better microbiota R2")
print(f"  3. Use the optimized data for downstream analysis")
print(f"  4. Generate plots and visualizations from the trained model")

def check_memory_usage():
    """Check current memory usage"""
    try:
        import psutil
        memory = psutil.virtual_memory()
        return f"Memory usage: {memory.percent:.1f}% ({memory.used/1e9:.1f}GB used / {memory.total/1e9:.1f}GB total)"
    except ImportError:
        return "Memory monitoring unavailable"

def reduce_features_by_variance(df, max_features=None, variance_threshold=0.01):
    """Reduce features by removing low-variance ones"""
    print(f"  Original shape: {df.shape}")
    
    # Separate sample column from data
    sample_col = df.columns[0]  # Assume first column is sample ID
    data_cols = df.columns[1:]
    
    # Calculate variance for numeric columns only
    numeric_data = df[data_cols].select_dtypes(include=[np.number])
    if len(numeric_data.columns) > 0:
        variances = numeric_data.var()
        high_var_cols = variances[variances > variance_threshold].index
        
        # Keep sample column and high-variance features
        keep_cols = [sample_col] + list(high_var_cols)
        df_filtered = df[keep_cols]
        print(f"  After variance filtering (>{variance_threshold}): {df_filtered.shape}")
    else:
        df_filtered = df
    
    # Further reduce if still too many features
    if max_features and df_filtered.shape[1] > max_features:
        # Keep sample column and top variance features
        data_cols = df_filtered.columns[1:]
        if len(data_cols) > max_features - 1:
            numeric_data = df_filtered[data_cols].select_dtypes(include=[np.number])
            if len(numeric_data.columns) > 0:
                variances = numeric_data.var().sort_values(ascending=False)
                top_features = variances.head(max_features - 1).index
                df_filtered = df_filtered[[sample_col] + list(top_features)]
                print(f"  After feature selection (top {max_features-1}): {df_filtered.shape}")
    
    return df_filtered

def optimize_data_types(df):
    """Convert to more memory-efficient data types"""
    df_opt = df.copy()
    
    # Convert float64 to float32
    float_cols = df_opt.select_dtypes(include=['float64']).columns
    df_opt[float_cols] = df_opt[float_cols].astype('float32')
    
    # Convert int64 to smaller int types where possible
    int_cols = df_opt.select_dtypes(include=['int64']).columns
    for col in int_cols:
        col_max = df_opt[col].max()
        col_min = df_opt[col].min()
        
        if col_min >= 0:  # Unsigned integers
            if col_max < 255:
                df_opt[col] = df_opt[col].astype('uint8')
            elif col_max < 65535:
                df_opt[col] = df_opt[col].astype('uint16')
            else:
                df_opt[col] = df_opt[col].astype('uint32')
        else:  # Signed integers
            if col_min > -128 and col_max < 127:
                df_opt[col] = df_opt[col].astype('int8')
            elif col_min > -32768 and col_max < 32767:
                df_opt[col] = df_opt[col].astype('int16')
            else:
                df_opt[col] = df_opt[col].astype('int32')
    
    return df_opt 
#                     "microbiota": microbiota_likelihood, 
#                     "metabolomics": "Normal"}
#     ),
#     mfl.TrainingOptions(
#         batch_size=43,  # Use all samples (small dataset)
#         max_epochs=1000  # Reduced from default
#     ),
# )reWarning)

# RESOURCE-OPTIMIZED DATA LOADING
print("=== MEMORY-OPTIMIZED DATA LOADING ===")
print(check_memory_usage())

print("\nLoading and optimizing data from CSV files...")

# Load the three datasets with memory optimization
lipidomics_df = pd.read_csv("Data/lipidomics_data_case.csv")
print(f"Loaded lipidomics: {lipidomics_df.shape}")

microbiota_df = pd.read_csv("Data/microbiota_SOFA_case_tumoral.csv")
print(f"Loaded microbiota: {microbiota_df.shape}")

metabolomics_df = pd.read_csv("Data/Metabolomics_data_case.csv")
print(f"Loaded metabolomics: {metabolomics_df.shape}")

print(f"\nOriginal dataset overview:")
print(f"Lipidomics data shape: {lipidomics_df.shape}")
print(f"Microbiota data shape: {microbiota_df.shape}")
print(f"Metabolomics data shape: {metabolomics_df.shape}")

# FEATURE REDUCTION FOR MEMORY EFFICIENCY
print("\n=== REDUCING FEATURES FOR MEMORY EFFICIENCY ===")


print("Reducing lipidomics features...")
lipidomics_df = reduce_features_by_variance(lipidomics_df, max_features=1000, variance_threshold=0.0)
lipidomics_df = optimize_data_types(lipidomics_df)

# Microbiota robustness patch: zero-filtering, less aggressive feature reduction, integer non-negative counts
print("Robustly reducing microbiota features...")
sample_col = microbiota_df.columns[0]

# Remove all-zero and near-zero features (columns)
feature_nonzero = (microbiota_df.iloc[:,1:] != 0).sum(axis=0)
feature_threshold = max(1, int(0.01 * microbiota_df.shape[0]))  # keep features with >1% nonzero
keep_features = feature_nonzero[feature_nonzero > feature_threshold].index
microbiota_df = pd.concat([microbiota_df[[sample_col]], microbiota_df.loc[:, keep_features]], axis=1)

# Remove all-zero and near-zero samples (rows)
sample_nonzero = (microbiota_df.iloc[:,1:] != 0).sum(axis=1)
sample_threshold = max(1, int(0.01 * (microbiota_df.shape[1]-1)))  # >1% nonzero features
rows_nz = sample_nonzero > sample_threshold
microbiota_df = microbiota_df.loc[rows_nz].reset_index(drop=True)

# Less aggressive feature reduction (keep all after cleaning)
microbiota_df = reduce_features_by_variance(microbiota_df, max_features=1000, variance_threshold=0.0)

# Ensure integer non-negative counts

# Only clip numeric columns (all except sample_col)
for col in microbiota_df.columns[1:]:
    microbiota_df[col] = pd.to_numeric(microbiota_df[col], errors='coerce').fillna(0)
    microbiota_df[col] = microbiota_df[col].clip(lower=0)
microbiota_df = optimize_data_types(microbiota_df)

# Print summary statistics for microbiota data
print("\n=== MICROBIOTA DATA SUMMARY BEFORE MODELING ===")
microbiota_numeric = microbiota_df.iloc[:,1:]
print(f"Shape: {microbiota_numeric.shape}")
print(f"Mean (per feature): {microbiota_numeric.mean().mean():.2f}")
print(f"Variance (per feature): {microbiota_numeric.var().mean():.2f}")
print(f"% zeros: {100 * (microbiota_numeric == 0).sum().sum() / microbiota_numeric.size:.2f}%")
print(f"Min: {microbiota_numeric.min().min()}, Max: {microbiota_numeric.max().max()}")

# Check for batch effects (simple: print mean per sample)
print("Mean per sample (first 10):", microbiota_numeric.mean(axis=1).head(10).tolist())

print("Reducing metabolomics features...")
metabolomics_df = reduce_features_by_variance(metabolomics_df, max_features=60, variance_threshold=0.001)
metabolomics_df = optimize_data_types(metabolomics_df)

print(f"\nOptimized dataset overview:")
print(f"Lipidomics data shape: {lipidomics_df.shape}")
print(f"Microbiota data shape: {microbiota_df.shape}")
print(f"Metabolomics data shape: {metabolomics_df.shape}")

# Force garbage collection
gc.collect()
print(f"\nAfter optimization: {check_memory_usage()}")

# Display first few rows and columns for each dataset
print("\n=== LIPIDOMICS DATA ===")
print("Columns:", lipidomics_df.columns.tolist()[:10])  # First 10 columns
print("First 5 rows:")
print(lipidomics_df.head())

print("\n=== MICROBIOTA DATA ===")
print("Columns:", microbiota_df.columns.tolist()[:10])  # First 10 columns
print("First 5 rows:")
print(microbiota_df.head())

print("\n=== METABOLOMICS DATA ===")
print("Columns:", metabolomics_df.columns.tolist()[:10])  # First 10 columns
print("First 5 rows:")
print(metabolomics_df.head())

# Check for common sample identifiers
print("\n=== SAMPLE INFORMATION ===")
# Assuming first column contains sample IDs (adjust if needed)
lipidomics_samples = set(lipidomics_df.iloc[:, 0])
microbiota_samples = set(microbiota_df.iloc[:, 0])
metabolomics_samples = set(metabolomics_df.iloc[:, 0])

print(f"Lipidomics samples: {len(lipidomics_samples)}")
print(f"Microbiota samples: {len(microbiota_samples)}")
print(f"Metabolomics samples: {len(metabolomics_samples)}")

# Find common samples across all datasets
common_samples = lipidomics_samples & microbiota_samples & metabolomics_samples
print(f"Common samples across all datasets: {len(common_samples)}")

# Check data types and missing values
print("\n=== DATA QUALITY CHECK ===")
for name, df in [("Lipidomics", lipidomics_df), ("Microbiota", microbiota_df), ("Metabolomics", metabolomics_df)]:
    print(f"\n{name}:")
    print(f"  Data types: {df.dtypes.value_counts().to_dict()}")
    print(f"  Missing values: {df.isnull().sum().sum()}")
    print(f"  Numeric columns: {len(df.select_dtypes(include=[np.number]).columns)}")

# Create MuData object for MOFAFLEX
print("\n=== PREPARING DATA FOR MOFAFLEX ===")

# Fix metabolomics column names issue
print("Fixing metabolomics data format...")
print("Metabolomics columns:", metabolomics_df.columns.tolist()[:5])
print("Metabolomics index:", metabolomics_df.index.tolist()[:5])

if 'Sample",Pyruvic acid"' in metabolomics_df.columns:
    # The first column is malformed - it contains sample IDs as the index
    # The sample IDs are actually in the index (M04, M05, etc.)
    metabolomics_df = metabolomics_df.reset_index()
    metabolomics_df = metabolomics_df.rename(columns={'index': 'Sample'})
    print("Fixed metabolomics sample IDs from index")
else:
    # Check if Sample column exists, if not create it from index
    print("No Sample column found, checking if sample IDs are in index...")
    if any(str(idx).startswith('M') for idx in metabolomics_df.index):
        metabolomics_df = metabolomics_df.reset_index()
        metabolomics_df = metabolomics_df.rename(columns={'index': 'Sample'})
        print("Created Sample column from index")
    else:
        print("Creating Sample column manually...")
        metabolomics_df.insert(0, 'Sample', [f'M{i+4:02d}' for i in range(len(metabolomics_df))])

print("After fixing - Columns:", metabolomics_df.columns.tolist()[:3])

# Update sample sets after fixing
lipidomics_samples = set(lipidomics_df['Sample'].astype(str))
microbiota_samples = set(microbiota_df['Sample'].astype(str))
metabolomics_samples = set(metabolomics_df['Sample'].astype(str))

print("Sample ID formats:")
print(f"Lipidomics: {sorted(list(lipidomics_samples))[:5]}")
print(f"Microbiota: {sorted(list(microbiota_samples))[:5]}")
print(f"Metabolomics: {sorted(list(metabolomics_samples))[:5]}")

# Find common samples across all datasets
common_samples = lipidomics_samples & microbiota_samples & metabolomics_samples
print(f"Common samples after fixing: {len(common_samples)}")

# Create AnnData objects for each modality
# Set sample column as index for each dataframe

lipidomics_data = lipidomics_df.set_index('Sample')
microbiota_data = microbiota_df.set_index('Sample')
metabolomics_data = metabolomics_df.set_index('Sample')

# Microbiota: ensure non-negative integer counts

# Log1p transform to reduce skew and prevent large magnitude issues
import numpy as np
from sklearn.preprocessing import StandardScaler
microbiota_data_log = np.log1p(microbiota_data)

# Optionally: standard scale across features
scaled_microbiota = pd.DataFrame(
    StandardScaler().fit_transform(microbiota_data_log),
    index=microbiota_data_log.index,
    columns=microbiota_data_log.columns
)

# Use this as input to AnnData
microbiota_adata = ad.AnnData(
    X=scaled_microbiota.values.astype(np.float32),
    obs=pd.DataFrame(index=scaled_microbiota.index),
    var=pd.DataFrame(index=scaled_microbiota.columns)
)

# Check if microbiota data is overdispersed

print("\n=== CHECKING MICROBIOTA DATA DISTRIBUTION ===")

mean_of_means = microbiota_data.mean(axis=0).mean()
mean_of_vars = microbiota_data.var(axis=0).mean()
vmr = mean_of_vars / max(mean_of_means, 1e-8)
print(f"Mean of means: {mean_of_means:.2f}, Mean of variances: {mean_of_vars:.2f}, VMR: {vmr:.2f}")
if vmr > 1.5:
    print("Data appears overdispersed (VMR > 1.5) - Negative Binomial likelihood recommended")
    microbiota_likelihood = "NegativeBinomial"
else:
    print("Data not heavily overdispersed (VMR <= 1.5) - Poisson likelihood appropriate")
    microbiota_likelihood = "Poisson"

print("\n=== DATA PREPARATION COMPLETED ===")
print(f"✅ Successfully processed 3 datasets with {len(common_samples)} matching samples")
print(f"✅ Recommended likelihood for microbiota: {microbiota_likelihood}")
print("✅ Data is ready for MOFAFLEX analysis!")

# Data summary
print(f"\nFinal data shapes:")
print(f"   - Lipidomics: {lipidomics_data.shape[0]} samples × {lipidomics_data.shape[1]} features")
print(f"   - Microbiota: {microbiota_data.shape[0]} samples × {microbiota_data.shape[1]} features") 
print(f"   - Metabolomics: {metabolomics_data.shape[0]} samples × {metabolomics_data.shape[1]} features")

# Save processed data
print("\nSaving processed data...")
lipidomics_data.to_csv("processed_lipidomics_data.csv")
microbiota_data.to_csv("processed_microbiota_data.csv") 
metabolomics_data.to_csv("processed_metabolomics_data.csv")
print("✅ Processed data saved to CSV files")

print(f"\n🎯 MOFAFLEX Configuration Recommendations:")
print(f"   - Lipidomics: Normal likelihood (continuous data)")
print(f"   - Microbiota: {microbiota_likelihood} likelihood (count data)")
print(f"   - Metabolomics: Normal likelihood (continuous data)")
print(f"   - Suggested factors: 5-10 (start small)")
print(f"   - Batch size: 43 (full dataset)")

print("\n=== CREATING MOFAFLEX MODEL ===")

# Create AnnData objects with memory-optimized data types
print("Creating memory-optimized AnnData objects...")

lipidomics_adata = ad.AnnData(
    X=lipidomics_data.values.astype(np.float32),  # Use float32 instead of float64
    obs=pd.DataFrame(index=lipidomics_data.index),
    var=pd.DataFrame(index=lipidomics_data.columns)
)

microbiota_adata = ad.AnnData(
    X=microbiota_data.values.astype(np.float32),  # Use float32 for MOFAFLEX compatibility
    obs=pd.DataFrame(index=microbiota_data.index),
    var=pd.DataFrame(index=microbiota_data.columns)
)

metabolomics_adata = ad.AnnData(
    X=metabolomics_data.values.astype(np.float32),  # Use float32 instead of float64
    obs=pd.DataFrame(index=metabolomics_data.index),
    var=pd.DataFrame(index=metabolomics_data.columns)
)

# --- Guiding variable integration ---
# Load guiding variable table from Data/ directory
tumor_location_df = pd.read_csv("Data/Metadata_Tumor_Location.csv")
    # Pad sample IDs: M4 → M04, M9 → M09
tumor_location_df["Sample"] = tumor_location_df["Sample"].apply(lambda x: f"M{int(x[1:]):02}")
tumor_location_df.set_index("Sample", inplace=True)

# Convert one-hot columns to single categorical label
# (assumes columns are one-hot encoded for each location)
tumor_location_series = tumor_location_df.idxmax(axis=1).str.replace("Tumor_", "")

# Add guiding variable to AnnData objects

# Fill NaNs and cast to category for tumor_location
for adata in [lipidomics_adata, microbiota_adata, metabolomics_adata]:
    adata.obs["tumor_location"] = tumor_location_series.reindex(adata.obs.index).fillna("Unknown").astype("category")
    print(f"Number of NaNs in tumor_location for {adata}:", adata.obs["tumor_location"].isna().sum())
    print(f"Categories: {adata.obs['tumor_location'].unique()}")

# --- Validation and diagnostics ---
meta_samples = set(tumor_location_series.index)
for modality, adata in [
    ("lipidomics", lipidomics_adata),
    ("microbiota", microbiota_adata),
    ("metabolomics", metabolomics_adata)
]:
    missing = set(adata.obs.index) - meta_samples
    if missing:
        print(f"[WARNING] {modality}: {len(missing)} samples missing tumor metadata: {sorted(list(missing))[:10]}")
    print(f"{modality}: tumor_location value counts:")
    try:
        print(adata.obs['tumor_location'].value_counts(dropna=False), "\n")
    except Exception as _e:
        print(f"Could not print tumor_location counts for {modality}: {_e}\n")

# Quick integrity checks for metadata one-hot format
try:
    one_hot_df = tumor_location_df.copy()
    if 'Sample' in one_hot_df.columns:
        one_hot_df = one_hot_df.drop(columns=['Sample'], errors='ignore')
    row_sums = one_hot_df.sum(axis=1)
    if (row_sums != 1).any():
        bad = (row_sums != 1).sum()
        print(f"[WARNING] Metadata one-hot integrity: {bad} rows do not sum to 1 (expected one-hot). Check tumor_location_df.sum(axis=1).head() for examples.")
except Exception:
    # If tumor_location_df is not one-hot (e.g., already a single column), skip this check
    pass

# --- DataOptions for MOFAFLEX ---
from mofaflex import DataOptions

data_opts = DataOptions(
    guiding_vars_obs_keys="tumor_location",
    subset_var=None,
    remove_constant_features=True,
    scale_per_group=False,
    plot_data_overview=False
)

# Create MuData object

# === Microbiota constant/all-zero feature/sample checks and filtering ===
from scipy.sparse import issparse
X_micro = microbiota_adata.X.toarray() if issparse(microbiota_adata.X) else microbiota_adata.X
feature_var = np.var(X_micro, axis=0)
zero_var_features = np.sum(feature_var == 0)
print(f"Number of constant features in microbiota: {zero_var_features}")
row_sums = X_micro.sum(axis=1)
zero_rows = np.sum(row_sums == 0)
print(f"Number of all-zero rows in microbiota: {zero_rows}")
col_sums = X_micro.sum(axis=0)
zero_cols = np.sum(col_sums == 0)
print(f"Number of all-zero columns in microbiota: {zero_cols}")
# Filter all-zero rows/columns
keep_rows = row_sums != 0
keep_cols = col_sums != 0
microbiota_adata = microbiota_adata[keep_rows, keep_cols].copy()

# Rebuild MuData with filtered microbiota_adata
mdata = md.MuData({
    'lipidomics': lipidomics_adata,
    'microbiota': microbiota_adata,
    'metabolomics': metabolomics_adata
})

print(f"Created optimized MuData object with {mdata.n_obs} observations and {len(mdata.mod)} modalities")
gc.collect()
print(f"Before model creation: {check_memory_usage()}")

# Create RESOURCE-OPTIMIZED MOFAFLEX model
print("\n=== CREATING RESOURCE-OPTIMIZED MOFAFLEX MODEL ===")

# Calculate optimal batch size (smaller for limited memory)
optimal_batch_size = len(common_samples)
print(f"  - Using optimal batch size: {optimal_batch_size} (full-batch for max memory usage)")

# Create MOFAFLEX model with user-specified options

# === Define model options ===

# Increase model capacity for microbiota
model_opts = mfl.ModelOptions(
    n_factors=20,
    likelihoods={
        "lipidomics": "Normal",
        "microbiota": "Normal",
        "metabolomics": "Normal"
    }
)

# === Define training options ===

training_opts = mfl.TrainingOptions(
    batch_size=optimal_batch_size,
    max_epochs=2000,
    seed=42
)

# === Create and initialize the model ===

# Print actual likelihood used for microbiota
print(f"Microbiota likelihood: {model_opts.likelihoods['microbiota']}")

model = mfl.MOFAFLEX(
    mdata,
    model_opts,
    training_opts,
    data_opts
)

# ============================
# 1️⃣ TRAINING PERFORMANCE
# ============================
print("\n=== INSPECTING TRAINING PERFORMANCE ===")
fig = mfl.pl.training_curve(model)
save_plot(fig, "training_curve.png", dpi=300)

# ============================
# 2️⃣ MODEL STRUCTURE INSPECTION
# ============================
print("\n=== CHECKING FACTOR INDEPENDENCE ===")
fig = mfl.pl.factor_correlation(model, figsize=(8, 8))
save_plot(fig, "factor_correlation_heatmap.png", dpi=300)

print("\n=== VARIANCE EXPLAINED BY GROUP AND VIEW ===")
fig1 = mfl.pl.variance_explained(model, group_by="group")
save_plot(fig1, "variance_explained_by_group.png", dpi=300)

fig2 = mfl.pl.variance_explained(model, group_by="view", figsize=(8, 5))
save_plot(fig2, "variance_explained_by_view.png", dpi=300)

# ============================
# 3️⃣ FACTOR EXPLORATION
# ============================
print("\n=== VISUALIZING LATENT FACTORS ===")
# Color factors by tumor location (guiding variable)
# Many mofaflex versions provide a scatter function to plot factors colored by covariates
try:
    # plot factor 1 vs 2 colored by tumor_location
    # configure desired plotting parameters (adjust these as needed)
    x_factor = 1
    y_factor = 2
    groups = None            # show all groups by default; set to a string or list to subset
    color_by = 'tumor_location'  # guiding covariate to color points
    shape_by = None          # set to 'tumor_location' to map shapes as well (may clutter)
    size = 2
    alpha = 0.9
    figsize = (6, 6)
    nrow = None
    ncol = None

    if hasattr(mfl.pl, 'factors_scatter'):
        # Call with the full parameter set supported by mofaflex.pl.factors_scatter
        try:
            fig = mfl.pl.factors_scatter(
                model,
                x=x_factor,
                y=y_factor,
                groups=groups,
                color=color_by,
                shape=shape_by,
                size=size,
                alpha=alpha,
                figsize=figsize,
                nrow=nrow,
                ncol=ncol,
            )
        except TypeError:
            # Some mofaflex versions may have a slightly different signature; try a minimal call
            fig = mfl.pl.factors_scatter(model, x_factor, y_factor, color=color_by, figsize=figsize)
    elif hasattr(mfl.pl, 'factors_covariate'):
        # older/newer API may provide a covariate plotting helper
        fig = mfl.pl.factors_covariate(model, covariate=color_by)
    else:
        # fallback: simple single-factor plot
        fig = mfl.pl.factor(model, factor=1)
except Exception:
    # As a last resort, try singular/plural names without coloring
    if hasattr(mfl.pl, 'factor'):
        fig = mfl.pl.factor(model, factor=1)
    elif hasattr(mfl.pl, 'factors'):
        fig = mfl.pl.factors(model)
    else:
        raise

save_plot(fig, "factors_colored_by_tumor_location.png", dpi=300)

# ============================
# 4️⃣ FEATURE LOADINGS (WEIGHTS)
# ============================
print("\n=== EXPLORING FEATURE WEIGHTS ===")

# (A) Full distribution of weights
mfl.pl.weights(model)
mfl.pl.weights(model, factors=1, n_features=15, pointsize=3)
mfl.pl.weights(model, factors=[1, 2, 3], figsize=(8, 6))

# Example: save weights for microbiota, Factor 1
fig = mfl.pl.weights(model, factors=1, views="microbiota")
save_plot(fig, "weights_factor1_microbiota.png", dpi=300)

# (B) Top-weighted features only
mfl.pl.top_weights(model, n_features=10)
mfl.pl.top_weights(model, factors=[1, 2], n_features=15)
mfl.pl.top_weights(model, views="microbiota", n_features=20, figsize=(6, 5))

# Example: save top features for lipidomics, Factor 1
fig = mfl.pl.top_weights(model, factors=1, views="lipidomics", n_features=15)
save_plot(fig, "top_weights_factor1_lipidomics.png", dpi=300)


# ============================
# 5️⃣ SUMMARY CHECKPOINTS
# ============================
print("\n=== SUMMARY CHECKS ===")
mfl.pl.training_curve(model)
mfl.pl.factor_correlation(model)
mfl.pl.variance_explained(model, group_by="view")
mfl.pl.weights(model, factors=1)
mfl.pl.top_weights(model, factors=1)

# Print R2 by view after training
try:
    r2 = model.get_r2()
    print("R2 by view:", r2)
except Exception as e:
    print("Could not compute R2 by view:", e)

# ============================
# 💾 SAVE TRAINED MODEL
# ============================
safe_model_save(model, "mofaflex_model_15factors.h5")
print("\n✅ Model and all plots successfully saved.")

# ============================
# 6️⃣ FACETED FACTOR PLOT (plotnine preferred, matplotlib fallback)
# ============================
try:
    # Get all factors as a single dataframe (posterior mean, mixed sparsity)
    factors_dict = model.get_factors(return_type='pandas', moment='mean', sparse_type='mix', ordered=True)

    # Assuming single group model (take first group)
    if isinstance(factors_dict, dict):
        Z = list(factors_dict.values())[0]
    else:
        Z = factors_dict

    # Ensure DataFrame
    if not isinstance(Z, pd.DataFrame):
        Z = pd.DataFrame(Z)

    # Prepare data
    Z_reset = Z.reset_index()
    first_col = Z_reset.columns[0]
    if first_col != 'Sample':
        Z_reset = Z_reset.rename(columns={first_col: 'Sample'})

    Z_long = Z_reset.melt(id_vars='Sample', var_name='Factor', value_name='Value')

    # Optional: sort factors numerically if needed
    def _factor_key(s):
        try:
            return int(s.split()[-1])
        except Exception:
            return s
    categories = sorted(Z_long['Factor'].unique(), key=_factor_key)
    Z_long['Factor'] = pd.Categorical(Z_long['Factor'].astype(str), ordered=True, categories=categories)

    out_plot = f"factors_faceted.png"
    if HAVE_PLOTNINE:
        try:
            plot = (
                pn.ggplot(Z_long, pn.aes(x='Sample', y='Value')) +
                pn.geom_point(size=2, color='steelblue') +
                pn.facet_wrap('~ Factor', scales='free_y', ncol=4) +
                pn.theme_bw() +
                pn.theme(
                    axis_text_x=pn.element_text(rotation=90, size=6),
                    figure_size=(14, 8),
                    subplots_adjust={'wspace': 0.3, 'hspace': 0.4}
                ) +
                pn.labs(title=f"MOFA Factors per Sample", y="Factor Value", x="Sample")
            )
            plot.save(out_plot, dpi=300)
            print(f"Saved faceted factors plot with plotnine to: {out_plot}")
        except Exception as e_plotnine:
            print(f"Plotnine plot failed, will try matplotlib fallback: {e_plotnine}")
            raise
    else:
        # Matplotlib/seaborn fallback: small subplots per factor
        n_factors = Z.shape[1]
        ncols = 4
        nrows = int(np.ceil(n_factors / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3 * nrows), sharex=False)
        axes = np.atleast_2d(axes)
        axes_flat = axes.flatten()
        for i, factor in enumerate(Z.columns):
            ax = axes_flat[i]
            if HAVE_SEABORN:
                sns.scatterplot(x='Sample', y=factor, data=Z_reset, ax=ax, s=10)
            else:
                ax.plot(Z_reset['Sample'].values, Z_reset[factor].values, linestyle='none', marker='o', markersize=4)
            ax.set_title(str(factor))
            ax.tick_params(axis='x', rotation=90, labelsize=6)
        # Turn off unused axes
        for j in range(n_factors, len(axes_flat)):
            axes_flat[j].axis('off')
        plt.tight_layout()
        fig.savefig(out_plot, dpi=300)
        plt.close(fig)
        print(f"Saved faceted factors plot with matplotlib fallback to: {out_plot}")

except Exception as e:
    print(f"\u274C Could not create faceted factors plot: {e}")

# ============================
# EXPORT & INSPECT WEIGHTS
# ============================
print("\n=== VISUALIZING TOP WEIGHTS ===")
mfl.pl.top_weights(model)

print("\n=== EXTRACT NUMERICAL WEIGHTS (PANDAS) ===")
weights = model.get_weights(return_type='pandas')

print("\n=== SAVING WEIGHTS PER VIEW TO CSV ===")
for view, df in weights.items():
    outname = f"weights_{view}.csv"
    df.to_csv(outname)
    print(f"Saved weights for {view} to {outname}")

print("\n=== TOP LIPIDS FOR FACTOR1 ===")
if 'lipidomics' in weights and len(weights['lipidomics'].columns) > 0:
    first_factor = weights['lipidomics'].columns[0]
    print(f"Top lipids for {first_factor}:")
    top_lipids = weights['lipidomics'][first_factor].nlargest(15)
    print(top_lipids)
else:
    print("Could not find any factors in lipidomics weights output")


# ============================
# 7️⃣ EXTRACT & SAVE FACTORS (Z) AND WEIGHTS (W) + PLOTS
# ============================
try:
    import os
    os.makedirs('MOFA_outputs', exist_ok=True)

    # Extract factors Z
    try:
        factors_dict = model.get_factors(return_type='pandas', moment='mean', sparse_type='mix', ordered=True)
        Z = list(factors_dict.values())[0] if isinstance(factors_dict, dict) else factors_dict
    except Exception:
        # fallback: try model.Z or helper
        if hasattr(model, 'Z'):
            Z = pd.DataFrame(model.Z)
        else:
            Z = get_factors_mofaflex(model, return_type='pandas')

    # Extract weights W
    try:
        weights_dict = model.get_weights(return_type='pandas', moment='mean', sparse_type='mix', ordered=True)
    except Exception:
        weights_dict = model.get_weights(return_type='pandas') if hasattr(model, 'get_weights') else {}

    # Save Z and W
    Z.to_csv('MOFA_outputs/factors_Z.csv')
    for view, w in weights_dict.items():
        w.to_csv(f"MOFA_outputs/weights_W_{view}.csv")
    print("✅ Saved Z and W matrices to 'MOFA_outputs/'")

    # Plot heatmap of Z
    try:
        plt.figure(figsize=(12, 6))
        if HAVE_SEABORN:
            sns.heatmap(Z, cmap='vlag', center=0, cbar_kws={'label': 'Factor Value'})
        else:
            plt.imshow(Z.values, aspect='auto', cmap='vlag', vmin=-np.max(np.abs(Z.values)), vmax=np.max(np.abs(Z.values)))
            plt.colorbar(label='Factor Value')
            plt.yticks(np.arange(Z.shape[0]), Z.index)
            plt.xticks(np.arange(Z.shape[1]), Z.columns, rotation=90)
        plt.title('MOFAFLEX Factors (Z) per Sample')
        plt.xlabel('Factor')
        plt.ylabel('Sample')
        plt.tight_layout()
        plt.savefig('MOFA_outputs/factor_heatmap_Z.png', dpi=300)
        plt.close()
    except Exception as e:
        print(f"Could not plot factor heatmap: {e}")

    # Plot heatmaps of weights per view
    for view, W in weights_dict.items():
        try:
            plt.figure(figsize=(12, 8))
            if HAVE_SEABORN:
                sns.heatmap(W, cmap='vlag', center=0, cbar_kws={'label': 'Weight Value'})
            else:
                plt.imshow(W.values, aspect='auto', cmap='vlag', vmin=-np.max(np.abs(W.values)), vmax=np.max(np.abs(W.values)))
                plt.colorbar(label='Weight Value')
                plt.yticks(np.arange(W.shape[0]), W.index, fontsize=6)
                plt.xticks(np.arange(W.shape[1]), W.columns, rotation=90, fontsize=6)
            plt.title(f"MOFAFLEX Weights (W) for view: {view}")
            plt.xlabel('Factor')
            plt.ylabel('Feature')
            plt.tight_layout()
            plt.savefig(f"MOFA_outputs/weight_heatmap_W_{view}.png", dpi=300)
            plt.close()
        except Exception as e:
            print(f"Could not plot weight heatmap for {view}: {e}")

except Exception as e_all:
    print(f"Could not extract/save Z and W or plot heatmaps: {e_all}")

