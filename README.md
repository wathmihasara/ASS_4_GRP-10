# ASS_4_GRP-10
#PART 1 
#Question 1
#Downloading the gene from the following code.

```{r}
URL="https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
download.file(URL, destfile = "gene_expression.tsv")
```

url: Saves the file's URL for future downloads. This should be changed to the actual web address hosting your gene_expression.tsv file.
destfile: Indicates the name and location on your local system where the downloaded file will be saved. In this case, the file with the name gene_expression.tsv will be saved in the current working directory.

#read the file

```{r}
gene_data <- read.table("gene_expression.tsv", header = TRUE, sep = "\t")
```
the above code read.table(): This is the R function that reads tabular data files, TSV, CSV, etc., into a data frame. These function arguments provide options for how the data are to be read and interpreted.
"gene_expression.tsv": It gives the name of the file to be read. It assumes that the file in is the current working directory, for other location full path need to be mentioned: for example, "path/to/ gene_expression.tsv".
header = TRUE: It means the first line of the text file contains column names.
sep = "\t" : This informs the reader that the tab character (\t) is utilized as the delimiter between values, a convention in TSV files. local system where the downloaded file will be saved. In this case, the file with the name gene_expression.tsv will be saved in the current working directory.

#displaying the first six genes

```{r}
head(gene_data)
```
(gene_data): The first six rows of the dataset are displayed by head(gene_data) to provide a brief overview of its organization and content.

Question 2

#calculating the mean
```{r}
str(gene_data)
gene_data$mean_expression <- rowMeans(gene_data[, sapply(gene_data, is.numeric)], na.rm = TRUE)
head(gene_data)
```
This code inspects the structure of the RNA-seq dataset, computes the mean expression across numeric columns for each gene-counts coming from several samples. The outcome of the calculated mean was kept in a new column called mean_expression.
str(gene_data): Shows the structure of gene_data. It gives information such as the row number corresponding to genes, the column number representing samples, the types of data, and the column names.
Using the sapply() function with the function is.numeric() to be executed on each column in gene_data. The output is a logical vector indicating in which columns numeric data are stored (e.g. RNA-seq counts).
gene_data[, sapply(gene_data, is.numeric)]: Excludes the only numeric columns for mean calculation. 
rowMeans():  computes the mean of each row (gene) across the numeric columns (samples).
na.rm = TRUE:  ensures that any missing values (NA) are ignored in the calculation.
gene_data$mean_expression <- ... : It calculates the mean expression for each gene and stores it in a new column called mean_expression in the gene_data data frame.
head(gene_data) : Display the first six rows of the updated dataframe which should include the new column mean_expression.

Question 3

# listing the genes with highest mean expression
```{r}
top_genes <- head(gene_data[order(-gene_data$mean_expression), ], 10)
top_genes
```
order() : returns the indices which would sort a vector in ascending fashion.
The use of the- before gene_data$mean_expression sorts the numbers in a descending fashion - this means from highest to lowest mean expression.
gene_data[order(-gene_data$mean_expression), ] : This reordered the rows of gene_data according to the sorted mean expression values so genes with highest mean are at top.
head() : It selects the first 10 rows of the reordered data frame i.e., the 10 genes having highest mean expression.
top_genes <-. : These 10 selected genes are stored in the data frame top_genes for further inspection.
top_genes : It prints the top 10 genes along with their expression values and mean expression.



