### FAQs

**What information does the App contain?**
- The app allows interactive display of different features (such as cell_type, tissue, donor etc) on the UMAP plots.
- The app is organized by the gene names. User can select one or multiple genes from the drop-down menu on the top-left corner. 
- Expression for selected genes is summarized using UMAP, heatmap, boxplot and barplots.
- Additionally, users can score expression correlation between the selected genes. 
- For several studies, the differential gene expression for each cell type is provided in the downloadable tables.


**Comparison of studies**
- User can select any two studies to compare the UMAP plots together.
- For the selected genes, expression of first gene is also compared in the adjacent UMAP plots.
    
**UMAP, t-SNE plots**
- "Dimensionality reduction is a powerful tool for machine learning practitioners to visualize and understand large, high dimensional datasets"
- [Refer the link below for details](https://pair-code.github.io/understanding-umap/#:~:text=In%20the%20simplest%20sense%2C%20UMAP,as%20structurally%20similar%20as%20possible.&text=In%20order%20to%20construct%20the,a%20%22fuzzy%20simplicial%20complex%22) 
- [Additionally, refer this!](https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668)


**What is a typical time to display contents in the app?**
- A typical time to extract all tables and plots for a given gene is <5 seconds. 
- The app is capable of handling and displaying millions of cells

**Whom shall I contact for more details or if the App breaks down?**
- If you see a gray screen wrap while requesting certain displays, this suggests that the app broke down!!
- Please refresh the page and navigate to the same request.
- It is possible that the processes saturate the memory/java-scokets etc. leading to the black-out! 
- If the problem persists, you can contact any member of the data science team.

**I am a developer, can I re-use the codes?**
- [Application code repository](https://gitlab.bayer.com/pdd_lab/pddsinglecell_shinyapp)
- The codes use licensed display functions and hence have a copyright! Thus parts of the visualization functions can not be used/copied/developed further! 
- Please contact Dhawal Jain (dhawal.jain@bayer.com) for details.

**Data Access**
- All datasets displayed in the app are referenced! Please visit the study details section! 
- Physical location of the raw-datasets and SQLite databases: /gpfs01/home/glxaw/data/scRNASeq_datasets (Bayer HPC) 
- Refer the [Excel-Sheet](documents/study_details.xlsx) for more details.

