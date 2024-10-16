# ASS_4_GRP-10

# PART 1 

# Question 1

# Downloading the gene from the following code.

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

# Top ten genes with the highest mean Expression
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

# Standard deviation

#sd
```{r}
sd(south[,1])
sd(south[,2])
```

#sd
```{r}
sd(north[,1])
sd(north[,2])
```

sd(): This function computes the standard deviation of values in the column specified.
south[, 1]: This actually is the notation for accessing the first column of the dataframe south. Let's assume this column contains the measurements of the circumference of trees for one particular year, say 2005.
south[, 2]: Likewise, the following line will compute the standard deviation of the second column in the south data frame-which one can guess contains the tree circumference measurements for another year, say 2020.


# Boxplot

```{r}
boxplot(north[,1],north[,2],ylab="tree circumference",names = c("Circumf_2005_cm","Circumf_2020_cm"))
```

#box plot
```{r}
boxplot(south[,1],south[,2],ylab="tree circumference",names = c("Circumf_2005_cm","Circumf_2020_cm"))
```

boxplot(): This function generates the boxplots for both southwest and northeast data.


Question 9 

# Calculate the mean growth over the last 10 years at each site

```{r}
north2 <- subset(csv_columns,Site==("northeast"))
```

```{r}
north2 <- subset(csv_columns,Site==("northeast"))
north3 <- subset(north2[,3:6])
colMeans(north3)
```
The first line filters the csv_columns data to include only the rows where the Site is "northeast." This filtered data is saved in a new variable called north2, which contains information specific to the northeast site. The colMeans(north3) line then computes the average values for the selected columns in north3. Ultimately, the calculated averages will be displayed in their respective columns.

```{r}
south2 <- subset(csv_columns,Site==("southwest"))
south3 <- subset(south2[,3:6])
colMeans(south3)
```
The first line filters the csv_columns data to include only the rows where the Site is "southwest." This filtered data is stored in a new variable called south2, which contains information specific to the southwest site. The second line selects columns 3 to 6 from south2 and saves this subset in south3. Finally, the colMeans(south3) line calculates the average values for the selected columns in south3. The resulting averages will be displayed in their respective columns.


```{r}
str(csv_columns)
```
This command will show the entire data frame stored in the variable csv_columns and can inspect the data, including the columns that contain tree circumference measurements.

# mean

subset(): It generates a new data frame, north/south, containing only the rows whose value in the column Site is "northeast” or “Southeast."
select = c("Circumf_2005_cm", argument Circumf_2020_cm: This argument selects the column that will remain in the new data frame; here, it selects columns corresponding to girth measurements for 2005 and 2020.
colMeans(): A function that calculates the mean of numeric columns in north data frame.
This will output a named numeric vector showing the mean circumference for Circumf_2005_cm and Circumf_2020_cm.

Question 10

# Use the t.test to estimate the p-value that the 10 year growth is different at the two sites.

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




# Part 2

# Question 1

# Intsalling seqinr
```{r}
install.packages("seqinr")
```
This code is used to install the seqinr package.

```{r}
library(seqinr)
```
This code is used to load the seqinr package to the library.
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
This code is used to summarize the number of coding sequences for both organisms Escherichia coli and Deinococcus radiodurans. This code organizes the data into table format.

# Printing the table

```{r}
print(coding_sequences_table)
```
We used this command to print the contents of the data.

# Question 2

# Loading the library
```{r}
library(seqinr)
```
The code reads an uncompressed FASTA file named ecoli_cds.fa.gz, which contains coding sequences for Escherichia coli, using the read.fa function from the seqinr  package.

# Reading the FASTA files
```{r}
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
deinococcus_cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
```
This code is used to read and load the coding sequence data for both Escherichia coli and Deinococcu radiodurans Fasta files.

# Calculating the lengths of coding sequences
```{r}
ecoli_lengths <- sapply(ecoli_cds, length)
deinococcus_lengths <- sapply(deinococcus_cds, length)
```
These codes are used to calculate and store the lengths of all coding sequences for both organisms.

# Summing the lengths
```{r}
ecoli_total_length <- sum(ecoli_lengths)
```
```{r}
deinococcus_total_length <- sum(deinococcus_lengths)
```
These codes used to summarize the lengths of all coding sequence for both organisms

# Summary table
```{r}
coding_dna_table <- data.frame(
  Organism = c("E. coli", "Deinococcus radiodurans"),
  Total_Coding_DNA_Length = c(ecoli_total_length, deinococcus_total_length)
)
```
This code is used to summarize the number of coding sequences for both organisms Escherichia coli and Deinococcus radiodurans. This code organizes the data into table format

# Table with the data
```{r}
print(coding_dna_table)
```
We used this command to print the contents of the data.


# Differances 
```{r}
if (ecoli_total_length > deinococcus_total_length) {
  coding_dna_comparison <- "E. coli has more total coding DNA than Deinococcus radiodurans."
} else if (ecoli_total_length < deinococcus_total_length) {
  coding_dna_comparison <- "Deinococcus radiodurans has more total coding DNA than E. coli."
} else {
  coding_dna_comparison <- "Both organisms have the same total coding DNA length."
}
```
This code is used to determine and make a comparison of the total coding DNA lengths of the two organisms, allowing for the interpretation of which organism has more DNA sequences

# Description
```{r}
print(coding_dna_comparison)
```
This code is used to present the result of the earlier comparison, providing clear and informative information regarding which organism has more total sequence or their lengths are equal.

# Question 3

# Loading the libraries
```{r}
library(seqinr)
library(ggplot2)
```
The code reads an uncompressed FASTA file named ecoli_cds.fa.gz, which contains coding sequences for Escherichia coli, using the read.fa function from the seqinr  package.

# Merging the lengths to plot the box plot
```{r}
length_data <- data.frame(
  Length = c(ecoli_lengths, deinococcus_lengths),
  Organism = c(rep("E. coli", length(ecoli_lengths)), rep("Deinococcus radiodurans", length(deinococcus_lengths)))
)
```
This code creates data that merges the lengths of coding sequences for both organisms, making it suitable for generating box plots to visually compare the distribution of sequence lengths between E. coli and Deinococcus radiodurans.

# The box plot
```{r}
ggplot(length_data, aes(x = Organism, y = Length)) +
  geom_boxplot() +
  labs(title = "Length Analysis of Coding Sequences", x = "Species", y = "Sequence Length (nucleotides)") +
  theme_minimal()
```
This code is used to generate a box plot that visually compares the lengths of coding sequences between E. coli and Deinococcus radiodurans.

# Calculating the mean and median lengths
```{r}
ecoli_mean <- mean(ecoli_lengths)
ecoli_median <- median(ecoli_lengths)
deinococcus_mean <- mean(deinococcus_lengths)
deinococcus_median <- median(deinococcus_lengths)

```
This code calculates and stores both the mean and median lengths of coding sequences for E. coli and Deinococcus radiodurans, providing measures that summarize the distribution of sequence lengths for each organism.

# Summary of mean and median
```{r}
summary_length_table <- data.frame(
  Organism = c("E. coli", "Deinococcus radiodurans"),
  Mean_Length = c(ecoli_mean, deinococcus_mean),
  Median_Length = c(ecoli_median, deinococcus_median)
)
```
This code prepares summary data the mean and median lengths for the coding sequences of E. coli and Deinococcus radiodurans, facilitating straightforward comparison and analysis of their coding sequence characteristics.

# Table
```{r}
print(summary_length_table)
```
This code is used to print the summary table, allowing to easily view and compare the mean and median lengths of coding sequences for both organisms.

# Describing the differences of mean
```{r}
mean_diff <- if (ecoli_mean > deinococcus_mean) {
  "E. coli has a longer mean coding sequence length."
} else if (ecoli_mean < deinococcus_mean) {
  "Deinococcus radiodurans has a longer mean coding sequence length."
} else {
  "Both organisms have the same mean coding sequence length."
}
```
This code is used to generate a comparison statement about the mean coding sequence lengths of the two organisms, allowing for a clear interpretation of which organism has a longer mean length or if they are equal

# Describing the median difference 
```{r}
median_diff <- if (ecoli_median > deinococcus_median) {
  "E. coli has a longer median coding sequence length."
} else if (ecoli_median < deinococcus_median) {
  "Deinococcus radiodurans has a longer median coding sequence length."
} else {
  "Both organisms have the same median coding sequence length."
}
```
This code is used to generate a comparison statement regarding the median coding sequence lengths of the two organisms. It provides a clear interpretation of which organism has a longer median length or if they are equal.

# Data (final)
```{r}
print(mean_diff)
print(median_diff)
```

These codes are used to print the results of the comparisons for both the mean and median coding sequence lengths, allowing to easily view and interpret the findings regarding the differences between the two organisms

# Question 4

# Bar plot for ecoli

```{r}
install.packages("ggplot2")
```
This code is used to download the ggplot2 

# Loading the libraries

```{r}
library(seqinr)
library(ggplot2)

```
The above code was used to load the installed packages seqinr and ggplot2.

# Reading the fasta files
```{r}
ecoli_cds <- seqinr::read.fasta("ecoli_cds.fa")
deinococcus_cds <- seqinr::read.fasta("deinococcus_radiodurans_cds.fa")
```
This code is used to read and load the coding sequence data for both Escherichia coli and Deinococcu radiodurans Fasta files.

# Calculating the nucleotide frequency

```{r}
nucleotide_frequency <- function(cds) {
  dna <- unlist(cds)
  all_sequences <- paste(dna, collapse = "")
  base_frequency <- table(strsplit(all_sequences, split = "")[[1]])
  return(as.data.frame(base_frequency))
}

```
This code is used to calculate  the frequency of nucleotide bases in a set of coding sequences.


# Calculating the nucleotide frequencies

```{r}
ecoli_nucleotide_freq <- nucleotide_frequency(ecoli_cds)
deinococcus_nucleotide_freq <- nucleotide_frequency(deinococcus_cds)

```

These codes are used to calculate and store the nucleotide frequencies for the coding sequences of both organisms.

# Bar plot for E-coli
```{r}
library(ggplot2)
ggplot(ecoli_nucleotide_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Nucleotide Frequencies in E. coli", 
       x = "Nucleotide", 
       y = "Frequency") +
  theme_minimal()

```
This code creates a bar plot that visualizes the frequencies of different nucleotides in E. coli, making it easy to interpret the composition of its sequences.

# Bar plot for deinococcus

```{r}
library(ggplot2)

ggplot(deinococcus_nucleotide_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "yellow") +
  labs(title = "Nucleotide Frequencies in Deinococcus radiodurans", 
       x = "Nucleotide", 
       y = "Frequency") +
  theme_minimal()

```
This code creates a bar plot that visualizes the frequencies of different nucleotides in Deinococcus radiodurans, allowing for easy interpretation of its sequence composition.

# calculate amino acid frequencies

```{r}
cal_amino_acid_frequency <- function(cds) {
  # Translate coding sequences into protein sequences
  protein_sequences <- sapply(cds, function(seq) translate(seq))
  all_proteins <- paste(protein_sequences, collapse = "")
  # Calculate frequency of each amino acid
  aa_freq <- table(strsplit(all_proteins, NULL)[[1]])
  return(as.data.frame(aa_freq))
}
```
This code is used to calculate the frequencies of amino acids in the protein sequences translated from a set of coding sequences. This information is valuable for understanding the composition of proteins.


# Calculating the amino acid frequency

```{r}
ecoli_aa_freq <- cal_amino_acid_frequency(ecoli_cds)
deinococcus_aa_freq <- cal_amino_acid_frequency(deinococcus_cds)

```
These codes are used to calculate the frequencies of amino acids derived from the protein sequences translated from the coding sequences of both organisms.

# plot for e-coli
```{r}
ggplot(ecoli_aa_freq, aes(x = Var1, y = Freq)) +
   geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Amino Acid Frequencies in E. coli", x = "Amino Acid", y = "Frequency") +
  theme_minimal()
```
This code is used to create a bar plot that visualizes the frequencies of different amino acids in Escherichia coli. This visualization allows for easy interpretation of protein composition.

# calculate amino acid frequencies

```{r}
cal_amino_acid_frequency <- function(cds) {
  # Translate coding sequences into protein sequences
  protein_sequences <- sapply(cds, function(seq) translate(seq))
  all_proteins <- paste(protein_sequences, collapse = "")
  # Calculate frequency of each amino acid
  aa_freq <- table(strsplit(all_proteins, NULL)[[1]])
  return(as.data.frame(aa_freq))
}
```
This code is used to calculate the frequencies of amino acids in the protein sequences translated from a set of coding sequences. This information is valuable for understanding the composition of proteins.

# Calculating the aa frequency

```{r}
ecoli_aa_freq <- cal_amino_acid_frequency(ecoli_cds)
deinococcus_aa_freq <- cal_amino_acid_frequency(deinococcus_cds)
```
The code invokes the cal_amino_acid_frequency function for both E. coli and Deinococcus. It processes ecoli_cds, converting the coding sequences into protein sequences and tallying the frequencies of each amino acid, which are then stored in the variable ecoli_aa_freq. Likewise, it uses deinococcus_cds to analyze the amino acid frequencies for Deinococcus, saving the results in deinococcus_aa_freq

# plot for Deinococcus

```{r}
ggplot(deinococcus_aa_freq, aes(x = Var1, y = Freq)) +
   geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Amino Acid Frequencies in deinococcus", x = "Amino Acid", y = "Frequency") +
  theme_minimal()
```
This code is used to create a bar plot that visualizes the frequencies of different amino acids in Deinococcus radiodurans. This visualization allows for easy interpretation of protein composition.

# Checking the above values by printing them on the interface to check whether it is done correctly

```{r}

print("comparison of nucleotide frequencies:")
print(merge(ecoli_nucleotide_freq, deinococcus_nucleotide_freq, by = "Var1", suffixes = c("_ecoli", "_deinococcus")))
```
This code is used to print information and compare the nucleotide frequencies of E. coli and Deinococcus radiodurans, displaying the results in a clear format for analysis. It aids in understanding the similarities and differences in nucleotide composition between both organisms.

# Displaying the values

```{r}
print("Amino Acid Frequency Comparison:")
print(merge(ecoli_aa_freq, deinococcus_aa_freq, by = "Var1", suffixes = c("_ecoli", "_deinococcus")))

```
This code is used to print the information and compare the amino acid frequencies of E. coli and Deinococcus radiodurans, displaying the results in a clear format for analysis. This comparison aids in understanding the similarities and differences in protein composition between both organisms.

# Question 5

```{r}
library(seqinr)
```
Loading the installed seqinr package


# data will be stored in codon_list

```{r}
codon_usage_list <- list()
```
This code sets up a list named codon_usage_list, which will be used to store codon usage data later in the analysis.

# codon usage for all sequences
```{r}
codon_usage_list <- lapply(cds, function(seq) {
  uco(seq, as.data.frame = TRUE)
})

```
This code is used to calculate the codon usage for each sequence in the cds list and saves the results in codon_usage_list. And indicates that the output should be returned as a data frame, making it easier to work with and analyze.

# Calculating RSCU for all sequences

```{r}
rscu_list <- lapply(cds, function(seq) {
  uco(seq, index = "rscu", as.data.frame = TRUE)
})

```
This code is used to calculate the RSCU for all sequences in the cds list and stores the results in rscu_list. And indicates that the results should be returned in a data frame format for easier manipulation.

# Check the data frames

```{r}
head(codon_usage_list[[1]])
head(rscu_list[[1]])
```
This code is used to check the first few entries of the codon usage and RSCU data frames for the first sequence in their respective lists. This is useful for validating the calculations and ensuring the data is structured correctly.

# Checking the column names 
```{r}
print(names(codon_usage_list[[1]]))
print(names(rscu_list[[1]]))
```
This code is used to check and display the column names of the codon usage and RSCU data frames for the first sequence in their respective lists.


# Combine codon usage and RSCU
```{r}
combined_data_list <- mapply(function(cu, rscu) {
  merge(cu, rscu, by = "codon", suffixes = c("_usage", "_rscu"))
}, codon_usage_list, rscu_list, SIMPLIFY = FALSE)
```
This code is used to combine codon usage and RSCU data for all sequences by merging the corresponding data frames. The combined data is stored in, facilitating easier analysis and comparison of codon usage and RSCU metrics across different sequences.

# Combining all the results

```{r}
combined_codon_usage <- do.call(rbind, codon_usage_list)
```
This code is used to consolidate all individual codon usage data frames into one comprehensive data frame, combined_codon_usage, facilitating streamlined analysis and handling of the combined codon usage information


# summarise all the codon usage and RSCU for all the sequences
# Structure of the combined data frame

```{r}
head(combined_data_list[[1]])

```
This code is used to display the first few rows of the data frame stored in the first element of the sequence. This helps quickly inspect the structure and contents of the dataset.



# Summary of codon usage frequency
```{r}
summary_usage <- aggregate(freq_usage ~ codon, data = combined_data_list[[1]], sum)

```
This code computes the overall codon usage for every codon in the dataset. It organizes the data based on the codon column and totals the codon_usage values associated with each codon. The output, saved in the total_codon_usage variable, provides a summarized view of codon usage throughout the sequences.


# Summary of RSCU values
```{r}
summary_rscu <- aggregate(RSCU_usage ~ codon, data = combined_data_list[[1]], sum)

```
This code is used to calculate the total Relative Synonymous Codon Usage (RSCU) for each codon in the dataset. It groups the data by the codon column and sums the RSCU_usage values for each codon. The result, stored in the summary_rscu variable, provides a summarized view of RSCU values across the sequences.


# Summary of the above data
```{r}
print(summary_usage)
print(summary_rscu)
```
 These codes are used to print the contents of the summary usage and summary rscu data. This allows to review the total codon usage and the summed RSCU values for each codon, providing insights into how frequently each codon is used and its relative usage compared to expected values across the dataset.


# Visualization of the output of condon usage frequency
```{r}
ggplot(summary_usage, aes(x = codon, y = freq_usage)) +
  geom_bar(stat = "identity", fill = "brown") +
  labs(title = "Total Codon Frequency", x = "Codon", y = "Frequency") +
  theme_minimal()
```
This code is used to create a bar plot that visualizes the total frequency of each codon from the summary usage data. And to identify patterns and differences among codons.


# Visualization of the output of RSCU values
```{r}
ggplot(summary_rscu, aes(x = codon, y = RSCU_usage)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Total RSCU Values", x = "Codon", y = "RSCU") +
  theme_minimal()
```
This code is used to visually represent the output of RSCU values. and to identify which codons are used more or less frequently in comparison to expected usage.

# Question 6

# read protein sequences from FASTA files

```{r}
deinococcus_cds <- read.fasta("deinococcus_radiodurans_cds.fa", seqtype = "AA", as.string = TRUE)
ecoli_cds <- read.fasta("ecoli_cds.fa", seqtype = "AA", as.string = TRUE)
```
This code  used to read the protein sequences from the FASTA file for Deinococcus radiodurans and stores them in the variable deinococcus_cds. This function performs the same operation for the E. coli protein sequences, storing them in ecoli_cds.

# convert the protein sequences read from the FASTA files into character vectors
```{r}
deinococcus_proteins <- unlist(lapply(deinococcus_cds, as.character))
ecoli_proteins <- unlist(lapply(ecoli_cds, as.character))
```
This code is used to convert the protein sequences that were read from the FASTA files into character vectors for easier manipulation and analysis in R.

# extracts k-mers 
```{r}
get_kmers <- function(sequences, k) {
  kmers <- unlist(lapply(sequences, function(seq) {
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  return(kmers)
}
```
This code used to systematically extract all possible contiguous subsequences of length k from the provided protein sequences

# extracts k-mers of lengths 3, 4, and 5 from the protein sequences of Deinococcus radiodurans and E. coli.

```{r}
kmers_deinococcus <- unlist(lapply(3:5, function(k) get_kmers(deinococcus_proteins, k)))
kmers_ecoli <- unlist(lapply(3:5, function(k) get_kmers(ecoli_proteins, k)))
```
This code is used to gather all possible 3-mers, 4-mers, and 5-mers from the protein sequences of both organisms

# counting the frequency of each k-mer extracted from the protein sequences

```{r}
deinococcus_kmer_counts <- as.data.frame(table(kmers_deinococcus))
ecoli_kmer_counts <- as.data.frame(table(kmers_ecoli))
```
This code is used to generate and store the frequency of each k-mer, which is crucial for understanding the distribution and occurrence of specific k-mers within the protein sequences of both organisms


# Renaming the columns of the data frames

```{r}
colnames(deinococcus_kmer_counts) <- c("kmer", "freq_deinococcus")
colnames(ecoli_kmer_counts) <- c("kmer", "freq_ecoli")
```
This code is used to rename the columns of the data frame deinococcus_kmer_counts and ecoli_kmer_counts for clarity and ease of use.


# display the first few rows for the gene of interest Deinococcus 

```{r}
print(head(deinococcus_kmer_counts))

```
This code is used to display the first few rows of the deinococcus_kmer_counts data frame. this code helps to verify the contents and structure of the k-mer frequency data for Deinococcus radiodurans

# display the first few rows for e-coli

```{r}
print(head(ecoli_kmer_counts))
```
This code is used to display the first few rows of the ecoli_kmer_counts data frame. It helps to assess the contents and structure of the k-mer frequency data for E. coli,


# total frequencies of k-mers 

```{r}
total_deinococcus <- sum(deinococcus_kmer_counts$freq_deinococcus)
total_ecoli <- sum(ecoli_kmer_counts$freq_ecoli)
```
This code provides the total k-mer frequencies in both datasets, which can be useful for comparing k-mer distributions and understanding the overall abundance of k-mers in the protein sequences.

# normalizes the k-mer frequencies 

```{r}
deinococcus_kmer_counts$freq_deinococcus <- deinococcus_kmer_counts$freq_deinococcus / total_deinococcus
ecoli_kmer_counts$freq_ecoli <- ecoli_kmer_counts$freq_ecoli / total_ecoli
```
This code is used to normalize the frequencies for a better comparison of k-mer frequencies between the both organism by converting raw counts into proportions

# Combining the counts
```{r}
combined_counts <- merge(deinococcus_kmer_counts, ecoli_kmer_counts, by = "kmer", all = TRUE)
combined_counts[is.na(combined_counts)] <- 0
combined_counts$ratio <- combined_counts$freq_deinococcus / combined_counts$freq_ecoli
```
This code is used to combine the k-mer frequency counts from both Deinococcus radiodurans and E. coli into a single data frame and calculate the ratio of their frequencies.

# Sorting the combined counts
```{r}
sorted_counts <- combined_counts[order(-combined_counts$ratio), ]
```

This code is used to sort the combined k-mer counts based on the calculated ratio in descending order.


# Top 10 overrepresented k-mers
```{r}
top_overrepresented <- sorted_counts[1:10, ]
```

This code is used to extract the top 10 k-mers that are represented in Deinococcus radiodurans compared to E. coli.


# Plot for Top 10 Overrepresented k-mers
```{r}
ggplot(top_overrepresented, aes(x = reorder(kmer, ratio), y = ratio)) +
 geom_bar(stat = "identity", fill = "beige") +
  coord_flip() +
  labs(title = "Top 10 Overrepresented k-mers", x = "K-mers", y = "Ratio (Deinococcus / E. coli)")

```
This code used to create a bar plot that visualizes the top 10 epresented k-mers in Deinococcus radiodurans compared to E. coli. And sets the title of the plot and labels for the x and y axes, providing context for what the plot represents.

# Sorts the combined counts
```{r}
sorted_underrepresented <- combined_counts[order(combined_counts$ratio), ]
```
This code is used to sort the combined k-mer counts based on the ratio in ascending order. The sorted data frame is assigned to the variable sorted_underrepresented, which now contains the k-mers arranged by their frequency ratios from lowest to highest.

# Top 10 underrepresented k-mers
```{r}
top_underrepresented <- sorted_underrepresented[1:10, ]
```
This code is used to extract the top 10 k-mers that are underrepresented in Deinococcus radiodurans compared to E. coli. It helps to identify the k-mers that are significantly less prevalent in Deinococcus radiodurans, which can be important for understanding differences in protein composition, function, or evolutionary adaptations between the two organisms.


# Plot for Top 10 Underrepresented k-mers
```{r}
ggplot(top_underrepresented, aes(x = reorder(kmer, ratio), y = ratio)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(title = "Top 10 Underrepresented k-mers", x = "K-mers", y = "Ratio (Deinococcus / E. coli)")
```
This code is used to create a bar plot that visualizes the top 10 underrepresented k-mers in Deinococcus radiodurans compared to E. coli.


