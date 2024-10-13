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

# Read File 

```{r}
gene_data <- read.table("gene_expression.tsv", header = TRUE, sep = "\t")
```
the above code read.table(): This is the R function that reads tabular data files, TSV, CSV, etc., into a data frame. These function arguments provide options for how the data are to be read and interpreted.
"gene_expression.tsv": It gives the name of the file to be read. It assumes that the file in is the current working directory, for other location full path need to be mentioned: for example, "path/to/ gene_expression.tsv".
header = TRUE: It means the first line of the text file contains column names.
sep = "\t" : This informs the reader that the tab character (\t) is utilized as the delimiter between values, a convention in TSV files. local system where the downloaded file will be saved. In this case, the file with the name gene_expression.tsv will be saved in the current working directory.

# Displaying the first six genes

```{r}
head(gene_data)
```
(gene_data): The first six rows of the dataset are displayed by head(gene_data) to provide a brief overview of its organization and content.

Question 2

# Calculating the Mean
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

# Top ten genes with highest mean Expression
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

Question 4

# The number of genes with a mean <10
```{r}
num_genes_low_expression <- sum(gene_data$mean_expression < 10)
num_genes_low_expression
```
The code counts how many genes in the gene_data dataframe have a mean expression value below 10. It uses the sum() function to evaluate each entry in the mean_expression column, determining how many of these values are less than 10. The count of these low-expression genes is stored in the variable num_genes_low_expression, providing the total number of genes with low expression levels.

Question 5

# A histogram plot of the mean values
```{r}
install.packages("ggplot2")
library(ggplot2)
top_genes[,1:5]
top_genes[,c(1,5)]
Data_set<-top_genes[,c(1,5)]
```
ggplot(): Initializes the plot with Data set and maps mean expression to the x-axis.
top_genes[, 1:5]: Shows the first 5 columns of `top genes`
top_genes[, c(1, 5)]: Extracts only the first and fifth column
Data_set <- top_genes[, c(1, 5)]: Extracts Gene ID and mean expression column

# Histogram
```{r}
Data_set<-top_genes[,c(1,5)]
hist(Data_set$mean_expression,xlab="Mean Expression",ylab="Frequency",main="histogram plot of the mean
```
Data_set <- top_genes[, c(1, 5)]: This line creates a new data frame, Data_set, comprising only two columns, namely Gene_ID and mean_expression columns.
hist(): Plots histogram of the values in the mean_expression column
xlab: x-axis was named.
ylab: y-axis was named.
main: The plot was titled.
col: The bars are color-filled in light blue.
border: Bars are outlined in black borders.

Question 6

# CSV download
```{r}
URL="https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
download.file(URL, destfile = "gene_expression.csv")
```
This code is executed based on previously completed GitHub assessments provided.

# Column Names
```{r}
csv_columns <- read.csv("/home/s224747674/assesment 4_grp 10/ASS_4_GRP-10/gene_expression.csv")
colnames(csv_columns)
```
read.csv(): This method reads a CSV document from the route indicated and generates a data frame based on its contents. The path here should be the full directory path where the file is located on your local machine. The resulting data frame is stored in a variable called csv columns.
colnames(): This function extracts the names of columns in the data frame

Question 7

# Mean and standard deviation of three circumference at the start and end of the study
```{r}
csv_columns
```
This command will show the entire data frame stored in the variable csv_columns and can inspect the data, including the columns that contain tree circumference measurements.

# Mean
```{r}
north <- subset(csv_columns,Site==("northeast"),c("Circumf_2005_cm","Circumf_2020_cm"))
colMeans(north)
north
south <- subset(csv_columns,Site==("southwest"),c("Circumf_2005_cm","Circumf_2020_cm"))
colMeans(south)
south
```
subset(): It generates a new data frame, north/south, containing only the rows whose value in the column Site is "northeast” or “Southeast."
select = c("Circumf_2005_cm", argument Circumf_2020_cm: This argument selects the column that will remain in the new data frame; here, it selects columns corresponding to girth measurements for 2005 and 2020.
colMeans(): A function that calculates the mean of numeric columns in north data frame.
This will output a named numeric vector showing the mean circumference for Circumf_2005_cm and Circumf_2020_cm.

# Standard Deviation
```{r}
sd(north[,1])
sd(north[,2])
```
sd(): This function computes the standard deviation of values in the column specified.
south[, 1]: This actually is the notation for accessing the first column of the dataframe south. Let's assume this column contains the measurements of the circumference of trees for one particular year, say 2005.
south[, 2]: Likewise, the following line will compute the standard deviation of the second column in the south data frame-which one can guess contains the tree circumference measurements for another year, say 2020.

Question 8

# Mean and standard deviation of tree circumference at the start and end of the study
#mean 
```{r}
north <- subset(csv_columns,Site==("northeast"),c("Circumf_2005_cm","Circumf_2020_cm"))
colMeans(north)
north
```
#mean
```{r}
south <- subset(csv_columns,Site==("southwest"),c("Circumf_2005_cm","Circumf_2020_cm"))
colMeans(south)
south
```
This command will show the entire data frame stored in the variable csv_columns and can inspect the data, including the columns that contain tree circumference measurements.

#mean

subset(): It generates a new data frame, north/south, containing only the rows whose value in the column Site is "northeast” or “Southeast."
select = c("Circumf_2005_cm", argument Circumf_2020_cm: This argument selects the column that will remain in the new data frame; here, it selects columns corresponding to girth measurements for 2005 and 2020.
colMeans(): A function that calculates the mean of numeric columns in north data frame.
This will output a named numeric vector showing the mean circumference for Circumf_2005_cm and Circumf_2020_cm.


#Standard deviation

#sd
```{r}
sd(south[,1])
sd(south[,2])
```
#sd
```{r}
sd(north[,1])
sd(north[,2])

sd(): This function computes the standard deviation of values in the column specified.
south[, 1]: This actually is the notation for accessing the first column of the dataframe south. Let's assume this column contains the measurements of the circumference of trees for one particular year, say 2005.
south[, 2]: Likewise, the following line will compute the standard deviation of the second column in the south data frame-which one can guess contains the tree circumference measurements for another year, say 2020.


#Boxplot

```{r}
boxplot(north[,1],north[,2],ylab="tree circumference",names = c("Circumf_2005_cm","Circumf_2020_cm"))
```
#box plot
```{r}
boxplot(south[,1],south[,2],ylab="tree circumference",names = c("Circumf_2005_cm","Circumf_2020_cm"))
```

boxplot(): This function generates the boxplot.


Northeast
    Circumf_2010_cm      	Circumf_2015_cm     	Circumf_2020_cm
11.288                         	24.516                                   	54.228


Southwest
    Circumf_2010_cm      	Circumf_2015_cm     	Circumf_2020_cm
           10.106          	21.316          	45.596

Question 9 

#Calculate the mean growth over the last 10 years at each site

```{r}
north2 <- subset(csv_columns,Site==("northeast"))
```

```{r}
north2 <- subset(csv_columns,Site==("northeast"))
north3 <- subset(north2[,3:6])
colMeans(north3)
```
```{r}
south2 <- subset(csv_columns,Site==("southwest"))
south3 <- subset(south2[,3:6])
colMeans(south3)
```
```{r}
str(csv_columns)
```
This command will show the entire data frame stored in the variable csv_columns and can inspect the data, including the columns that contain tree circumference measurements.

#mean

subset(): It generates a new data frame, north/south, containing only the rows whose value in the column Site is "northeast” or “Southeast."
select = c("Circumf_2005_cm", argument Circumf_2020_cm: This argument selects the column that will remain in the new data frame; here, it selects columns corresponding to girth measurements for 2005 and 2020.
colMeans(): A function that calculates the mean of numeric columns in north data frame.
This will output a named numeric vector showing the mean circumference for Circumf_2005_cm and Circumf_2020_cm.

Qustion 10

#Use the t.test to estimate the p-value that the 10 year growth is different at the two sites.

```{r}
csv_columns$growth_10_year <- csv_columns$Circumf_2020_cm - csv_columns$Circumf_2010_cm
t_test_result <- t.test(growth_10_year ~ Site, data = csv_columns)
print(t_test_result)
```

csv_columns$growth_10_year <- csv_columns$Circumf_2020_cm - csv_columns$Circumf_2010_cm : This line creates a new column in the csv columns data frame called growth_10_year that calculates the growth in the circumference of trees, this done by subtracting the circumference measured in the year 2010 from that in 2020.
t.test(): This is a function of t-test that will show the comparison of mean values between two or more groups.
growth_10_year ~ Site : This model that you want to test the differences in growth_10_year between the groups defined by the Site variable.
data = csv_columns: This argument identifies the data frame containing the variables of interest.
print(t_test_result): This line outputs the results of the t-test, which includes statistics such as the t-value, degrees of freedom, p-value, and confidence interval for the mean difference.




#Part 2
Question 1

# Intsalling seqinr
```{r}
install.packages("seqinr")
```
```{r}
library(seqinr)
```
The code reads an uncompressed FASTA file named ecoli_cds.fa.gz, which contains coding sequences for Escherichia coli, using the read.fa function from the seqinr  package.

# Downloading the gene e-coli
```{r}
library("R.utils")
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile="ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
list.files()
```
Library (“R.utils”) this code loads the R.utils package, which provides ability to decompress files. Then it defines the URL where the compressed FASTA file is hosted.the command download.file(URL,destfile="ecoli_cds.fa.gz") downloads the file from the specified URL  and saved it with the name ecoli_cds.fa.gz. the command  gunzip("ecoli_cds.fa.gz") decompresses the downloaded gunzipped file. Finally, list.files() lists the files in the current working directory, confirming that the FASTA file is downloaded.

# Downloading the gene of interest
```{r}
library("R.utils")
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_0_collection/deinococcus_radiodurans_r1_gca_000008565/cds/Deinococcus_radiodurans_r1_gca_000008565.ASM856v1.cds.all.fa.gz"
download.file(URL,destfile="deinococcus_radiodurans_cds.fa.gz")
gunzip("deinococcus_radiodurans_cds.fa.gz")
list.files()
```
Library (“R.utils”) this code loads the R.utils package, which provides ability to decompress files. Then it defines the URL where the compressed FASTA file is hosted.the command download.file(URL,destfile="deinococcus_radiodurans_cds.fa.gz") downloads the file from the specified URL  and saved it with the name radiodurans_cds.fa.gz. the command  gunzip("deinococcus_radiodurans_cds.fa.gz") decompresses the downloaded gunzipped file. Finally, list.files() lists the files in the current working directory, confirming that the FASTA file is downloaded.

# Loading the sequence for both organisms
```{r}
cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
```
This code is used for load sequences for both Deinococcus radiodurans and Escherichia coli into R

# Counting the codes for ecoli
```{r}
ecoli_count <- length(ecoli_cds)
```
This code is used to determine how many sequences are stored in ecoli_cds.

# Counting the codes for deinococcus and Read the FASTA file for Deinococcus radiodurans
```{r}
deinococcus_cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
```
This code is used to read and load the coding sequence data for the Deinococcu radiodurans into R.

```{r}
deinococcus_count <- length(deinococcus_cds)
```
This code is used to determine how many sequences are stored in Deinococcu radiodurans. 

# Print the count for Deinococcus radiodurans
```{r}
print(paste("Number of coding sequences in Deinococcus radiodurans:", deinococcus_count)
```
The code is used to display the count of coding sequence for Deinococcus radiodurans and make sure to earlier used code work properly.

# Creating a summary table
```{r}
coding_sequences_table <- data.frame(
  Organism = c("E. coli", "Deinococcus radiodurans"),
  Number_of_Coding_Sequences = c(ecoli_count, deinococcus_count)
)
```
intsalling seqinr
```{r}
install.packages("seqinr")
```
```{r}
library(seqinr)
```
The code reads an uncompressed FASTA file named ecoli_cds.fa.gz, which contains coding sequences for Escherichia coli, using the read.fa function from the seqinr  package.

# Downloading the gene e-coli
```{r}
library("R.utils")
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile="ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
list.files()
```
Library (“R.utils”) this code loads the R.utils package, which provides ability to decompress files. Then it defines the URL where the compressed FASTA file is hosted.the command download.file(URL,destfile="ecoli_cds.fa.gz") downloads the file from the specified URL  and saved it with the name ecoli_cds.fa.gz. the command  gunzip("ecoli_cds.fa.gz") decompresses the downloaded gunzipped file. Finally, list.files() lists the files in the current working directory, confirming that the FASTA file is downloaded.

# Downloading the gene of interest
```{r}
library("R.utils")
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_0_collection/deinococcus_radiodurans_r1_gca_000008565/cds/Deinococcus_radiodurans_r1_gca_000008565.ASM856v1.cds.all.fa.gz"
download.file(URL,destfile="deinococcus_radiodurans_cds.fa.gz")
gunzip("deinococcus_radiodurans_cds.fa.gz")
list.files()
```
Library (“R.utils”) this code loads the R.utils package, which provides ability to decompress files. Then it defines the URL where the compressed FASTA file is hosted.the command download.file(URL,destfile="deinococcus_radiodurans_cds.fa.gz") downloads the file from the specified URL  and saved it with the name radiodurans_cds.fa.gz. the command  gunzip("deinococcus_radiodurans_cds.fa.gz") decompresses the downloaded gunzipped file. Finally, list.files() lists the files in the current working directory, confirming that the FASTA file is downloaded.

# Loading the sequence for both organisms
```{r}
cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
```
This code is used for load sequences for both Deinococcus radiodurans and Escherichia coli into R

# Counting the codes for ecoli
```{r}
ecoli_count <- length(ecoli_cds)
```
This code is used to determine how many sequences are stored in ecoli_cds.

# Counting the codes for deinococcus and Read the FASTA file for Deinococcus radiodurans
```{r}
deinococcus_cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
```
This code is used to read and load the coding sequence data for the Deinococcu radiodurans into R.

```{r}
deinococcus_count <- length(deinococcus_cds)
```
This code is used to determine how many sequences are stored in Deinococcu radiodurans. 

# Print the count for Deinococcus radiodurans
```{r}
print(paste("Number of coding sequences in Deinococcus radiodurans:", deinococcus_count)
```
The code is used to display the count of coding sequence for Deinococcus radiodurans and make sure to earlier used code work properly.

# Creating a summary table
```{r}
coding_sequences_table <- data.frame(
  Organism = c("E. coli", "Deinococcus radiodurans"),
  Number_of_Coding_Sequences = c(ecoli_count, deinococcus_count)
)
```
This code used to summarize the number of coding sequences for both organisms Escherichia coli and Deinococcus radiodurans. This code organizes the data into table format.

# Printing the table

```{r}
print(coding_sequences_table)
```
We used this command to print the contents of the data.





 



