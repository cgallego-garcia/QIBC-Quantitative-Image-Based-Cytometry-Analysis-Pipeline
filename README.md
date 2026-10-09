# QIBC Analysis Pipeline for ImageJ/Fiji

[![DOI](https://img.shields.io/badge/DOI-10.1126%2Fsciadv.adz6211-blue.svg)](https://doi.org/10.1126/sciadv.adz6211) Sara Martín-Vírgala et al. ,Slow RNAPII elongation enhances naive pluripotency rewiring while maintaining high replication fork speed.Sci. Adv.12,eadz6211(2026).

This ImageJ/Fiji macro automates the Quantitative Image-Based Cytometry (QIBC) workflow. It is designed to intelligently process both 2D images and 3D Z-stacks. The pipeline performs nuclei segmentation via Cellpose (using the DAPI channel), applies morphological filtering to remove artifacts, and classifies cells into four populations (Double Positive, Double Negative, EdU+, and H2Ax+) based on user-defined intensity thresholds.

## 🛠️ Prerequisites and Dependencies

To run this macro, you need [Fiji (ImageJ)](https://imagej.net/software/fiji/) installed and the following components configured:

1. **Fiji Update Sites**:
   Open Fiji, go to `Help > Update...`, click on `Manage update sites`, and activate:
   * **IJPB-plugins** (provides MorphoLibJ for label processing and filtering).
   * **PTBIOP** (provides Cellpose wrapper integration in ImageJ).

2. **Cellpose Environment**:
   * You must have a local Conda/Mamba environment with `cellpose` installed.
   * The macro uses a pre-trained model (configured as `cpsam` by default, or your own custom model). *Note: Make sure to update the `model_path` variable in the code if you are using a specific custom model on your machine.*

## 🚀 Usage

1. Clone this repository or download the `.ijm` file.
2. Open Fiji, go to `Plugins > Macros > Run...`, and select the macro file.
3. A Graphical User Interface (GUI) will open where you must define:
   * **Input images folder**: The directory containing your images (`.nd2`, `.czi`, `.lsm`, `.tiff`).
   * **Cellpose environment path**: Path to your conda environment (e.g., `C:\Users\User\miniconda3\envs\cellpose`).
   * **Nuclei Area Selection**: Minimum and maximum nuclei sizes (in square pixels) to discard debris or clusters.
   * **Intensity levels**: Intensity thresholds to consider a cell positive for EdU and H2Ax.
   * **Background Subtraction**: "Rolling Ball" radius to subtract the background for each channel.

### How does it handle dimensions?
* **Single-plane images (2D):** Quantified directly after applying background subtraction.
* **Z-stacks (3D):** The macro automatically generates a **SUM projection** for robust quantitative measurements, and a **MAX projection** which is saved exclusively for qualitative visualization.

## 📁 Generated Outputs

The macro automatically creates a subfolder named `Resultados` inside your input directory containing:

* `Resultados_Contajes.xls`: A table with the absolute counts and percentages of the 4 populations for each image.
* `Resultados_Intensidad.xls`: Raw data table with the mean intensities of DAPI, EdU, and H2Ax for each individual cell.
* `*_ROIset_nuclei.zip`: A zip file containing the Regions of Interest (ROIs) for all segmented nuclei.
* `*_Classification_Map.tiff`: A flat image with the ROIs overlaid and color-coded according to their classification (Green: Double+, Red: Double-, Blue: EdU+, Orange: H2Ax+).
* `*_MAX_projection.tiff`: (Only generated if the input was a Z-stack).

## 🖋️ Authors
* **Carlos Gallego-Garcia** 
_Faculty of Experimental Sciences, Universidad Francisco de Vitoria, Madrid, Spain._
_Advanced Light Microscopy Facility (SMOA), Centro de Biología Molecular Severo Ochoa (CBM)._
