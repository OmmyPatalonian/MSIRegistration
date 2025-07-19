# MSI Histology Overlay Application

A sophisticated Shiny web application for visualizing mass spectrometry imaging (MSI) data with histology image overlays, specifically designed for DESI (Desorption Electrospray Ionization) mass spectrometry analysis.

## Features

### Core Functionality
- **MSI Data Visualization**: Multiple visualization modes including first ion, total ion current (TIC), and custom ion selection
- **Histology Overlay**: Real-time overlay of histology images on MSI data with adjustable transparency
- **Image Transformations**: Scale, rotate, and translate histology images for precise alignment
- **Mathematical Operations**: Perform complex ion calculations (ratios, differences, sums, etc.)

### Data Processing
- **File Format Support**: 
  - imzML/ibd files (standard MSI format)
  - RDS files (Cardinal-compatible objects)
- **Ion Selection**: Interactive data tables for selecting specific m/z values
- **Feature Summarization**: Automatic computation of mean intensities and frequency statistics

### Visualization Options
- **Color Schemes**: 20+ color palettes including Spectral, Viridis, Plasma, Inferno
- **Image Enhancement**: Contrast enhancement, smoothing, and normalization options
- **Spectrum Plotting**: Interactive mass spectra with PPM error calculations
- **Font Customization**: Scalable fonts and styling options

## Files

- `Histology_DESI_overlay (3).R` - Main Shiny application
- `plot_Card_overlay_NEW (3).R` - Plotting module with visualization functions
- `sample_histology.png` - Sample histology image for testing
- `sample_msi_data.rds` - Sample MSI dataset in Cardinal format

## Installation

### Prerequisites
```r
# Install required packages
install.packages(c("shiny", "png", "jpeg", "ggplot2", "grid", "abind", "tools", "DT"))

# Install Cardinal from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Cardinal")
```

### Running the Application
```r
# Clone the repository
git clone https://github.com/OmmyPatalonian/MSI.git
cd MSI

# Run the application
Rscript "Histology_DESI_overlay (3).R"
```

## Usage

1. **Upload Data**: 
   - Upload a histology image (PNG/JPEG)
   - Upload MSI dataset (.rds or .imzML + .ibd files)

2. **Align Images**:
   - Use transformation controls to align histology with MSI coordinate system
   - Adjust scale, rotation, and translation parameters

3. **Analyze Ions**:
   - Select visualization mode (First ion, TIC, or Custom)
   - Choose specific ions for analysis
   - Apply mathematical operations for complex calculations

4. **Customize Visualization**:
   - Adjust transparency slider for optimal overlay
   - Select color schemes and enhancement options
   - Configure plot dimensions and styling

5. **Save Settings**:
   - Download current transformation settings
   - Load previously saved configurations

## Technical Details

### Architecture
- **Modular Design**: Separate UI and server modules for maintainability
- **Reactive Programming**: Real-time updates based on user input
- **Error Handling**: Comprehensive error checking and user notifications

### Key Functions
- `plot_card_UI()`: User interface for MSI visualization controls
- `plot_card_server()`: Server logic for data processing and plotting
- `cpal()`: Color palette management with fallback options

### Dependencies
- **Cardinal**: Mass spectrometry imaging data processing
- **Shiny**: Web application framework
- **Grid/PNG/JPEG**: Image handling and overlay functionality
- **DT**: Interactive data tables

## Sample Data

The repository includes sample data for testing:
- **Sample MSI Data**: Simulated 20x20 pixel dataset with 20 peaks
- **Sample Histology**: Synthetic tissue-like image with circular structures

To generate fresh sample data, uncomment the `generate_sample_data()` call in the main application file.

## Troubleshooting

### Common Issues
1. **Large File Uploads**: Maximum file size is set to 500MB
2. **Memory Usage**: Large datasets may require increased R memory limits
3. **Package Installation**: Ensure all Bioconductor packages are properly installed

### Error Messages
- "No data available": Check file format and data integrity
- "Text to be written must be a length-one character vector": Usually resolved by restarting the application

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test with sample data
5. Submit a pull request

## License

This project is open source. Please cite appropriately if used in research.

## Contact

For questions or issues, please open a GitHub issue or contact the repository owner.

---

**Note**: This application is designed for research purposes and requires familiarity with mass spectrometry imaging concepts and Cardinal package usage.
