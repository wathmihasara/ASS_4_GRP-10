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
