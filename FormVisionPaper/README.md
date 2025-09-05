# FormVisionPaper

This directory contains scripts and notebooks for reproducing figures from:

**Seung, H.S.** *Predicting visual function by interpreting a neuronal wiring diagram.* Nature 634, 113–123 (2024). https://doi.org/10.1038/s41586-024-07953-5

## Files

### Main Scripts
- **`FiguresSpatial.jl`** - Generates spatial analysis figures (receptive fields, connectivity maps, projections)
- **`FiguresNonspatial.jl`** - Generates non-spatial figures (connectivity matrices, statistics)
- **`DataS2.jl`** - Supplementary data analysis

### Configuration
- **`config.jl`** - Configuration file defining TARGETDIR for output figures
- **`Project.toml`** - Julia project dependencies
- **`Manifest.toml`** - Dependency version lock file

## Usage

### As Julia Scripts
```bash
cd FormVisionPaper
julia --project=.
```

```julia
julia> include("FiguresSpatial.jl")  # Generate spatial figures
julia> include("FiguresNonspatial.jl")  # Generate non-spatial figures
```

### As Jupyter Notebooks
The main `.jl` files are configured with Jupytext, so they can be opened directly as Jupyter notebooks. First install and enable the Jupytext extension:

```bash
pip install jupytext
jupyter lab --generate-config  # if you don't have a config file
jupyter serverextension enable jupytext
```

Then you can open the scripts as notebooks:
```bash
jupyter notebook FiguresSpatial.jl
jupyter notebook FiguresNonspatial.jl
```

Jupytext will automatically convert between `.jl` and `.ipynb` formats while keeping the source files as Julia scripts in version control.

## Output Structure

Figures are saved to the `figures/` subdirectory organized by analysis type:

### Spatial Analysis Outputs
- **`figures/Dm3p/`** - Individual cell montages for Dm3p neurons
- **`figures/Dm3q/`** - Individual cell montages for Dm3q neurons  
- **`figures/Dm3v/`** - Individual cell montages for Dm3v neurons
- **`figures/TmY4/`** - Individual cell montages for TmY4 neurons
- **`figures/TmY9q/`** - Individual cell montages for TmY9q neurons
- **`figures/TmY9q⊥/`** - Individual cell montages for TmY9q⊥ neurons
- **`figures/LC10ev/`** - Individual cell montages for LC10ev neurons
- **`figures/LC15/`** - Individual cell montages for LC15 neurons
- **Main directory figures**: Summary plots (angles, aspects, projections, ERFs)

### Key Features

#### Spatial Analysis (`FiguresSpatial.jl`)
- **Connectivity Maps**: Tm1→Dm3 and Tm1→TmY receptive field analysis
- **ERF Analysis**: Extended receptive fields via disynaptic pathways
- **Individual Cell Montages**: Automated generation of pathway-specific receptive fields for each cell
- **Projection Analysis**: Hexagonal coordinate projections of connectivity patterns
- **LC Cell Analysis**: LC10ev and LC15 pathway analysis

#### Individual Cell Montage Generation
The script includes an optimized `generate_target_montages()` function that:
- Finds strongest monosynaptic and disynaptic pathways to each target cell type
- Generates receptive field maps for each pathway
- Creates individual PDF montages for every cell of the target type
- Organizes output into cell-type-specific subdirectories

## Data Dependencies

The scripts automatically download required FlyWire Codex data via DataDeps.jl on first run. Data is cached locally for subsequent analyses.

## Requirements

- Julia 1.10+ 
- OpticLobe.jl package
- Dependencies listed in Project.toml (Luxor, Plots, ProgressMeter, ColorSchemes, etc.)
- Jupytext extension (for notebook functionality)

## Notes

- Scripts use hexagonal coordinate systems for spatial analysis
- Side filtering restricts analysis to right hemisphere by default
- Error handling includes missing value propagation for incomplete connectivity data
- Figures match those published in the Nature paper