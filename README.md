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
