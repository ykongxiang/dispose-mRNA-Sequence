library(readxl)
library(tidyverse)
library(dplyr)
library(openxlsx)
oct_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_blastx.xlsx")
colnames(oct_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
oct_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx")
oct_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_prf.xlsx")
oct_mrna$DNA_seqid <- gsub("Contig.*", "", oct_mrna$DNA_seqid)
va_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_blastx.xlsx")
colnames(va_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
va_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_mrna.xlsx")
va_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_prf.xlsx")

va_prf_mrna <- merge(va_prf, va_mrna, by = "DNA_seqid", all = TRUE)
oct_prf_mrna <- merge(oct_prf, oct_mrna, by ='DNA_seqid', all = TRUE)
#删除contig标注 使DNA信息纯净

na_indices <- which(is.na(oct_prf_mrna$Sequence))
for (i in seq_along(na_indices)) {
  na_index <- na_indices[i]
  for (j in seq(na_index, nrow(oct_prf_mrna))) {
    if (!is.na(oct_prf_mrna$Sequence[j])) {
      oct_prf_mrna$Sequence[na_index] <- oct_prf_mrna$Sequence[j]
      break 
    }
  }
}
na_rows <- which(is.na(oct_prf_mrna$Strand))
oct_prf_mrna <- oct_prf_mrna %>%
  group_by(Sequence) %>%
  mutate(HasNonNA = any(!is.na(Strand))) %>%
  ungroup()
oct_prf_mrna <- oct_prf_mrna %>%
  filter(!(is.na(Strand) & oct_prf_mrna$HasNonNA))

write.xlsx(va_prf_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/va_prf_mrna.xlsx')
write.xlsx(oct_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx')



if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
library(Biostrings)



extract_sequences <- function(mrna_file) {
  mrna_data <- read_excel(mrna_file)
  mrna_data <- na.omit(mrna_data)  
  mrna_data$FS_start <- as.numeric(mrna_data$FS_start)
  mrna_data$FS_end <- as.numeric(mrna_data$FS_end)
  convert_to_positive_strand <- function(sequence, strand) {
    if (strand == '-') {
      return(reverseComplement(DNAString(sequence)))
    } else {
      return(sequence)
    }
  }
  
  extract_and_reverse_substring <- function(text, start, end, strand) {
    reverse_needed <- (strand == '-')
    substring <- substr(text, start-2, end + 1)
    if (reverse_needed) {
      text <- reverseComplement(DNAString(text))
      substring <- substr(text,)
        return(substring)
    
  }
  seqs <- sapply(1:nrow(mrna_data), function(i) extract_and_reverse_substring(
    text = mrna_data$Sequence[i],
    start = mrna_data$FS_start[i],
    end = mrna_data$FS_end[i],
    strand = mrna_data$Strand[i]
  ))
  return(seqs)
}


# 调用主函数
a <- extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx')
mrna_file <- '/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx'
mrna_data1 <- read_excel(mrna_file)
mrna_data1 <- na.omit(mrna_data1)



# 计算每个字符串的个数
string_counts <- table(a)

# 计算总数
total_count <- sum(string_counts)

# 计算每种字符串的频率
string_frequencies <- string_counts / total_count

data <- data.frame(
  String = names(string_frequencies),
  Frequency = as.numeric(string_frequencies)
)

# 对频率进行降序排序
data <- data[order(-data$Frequency), ]
first_six <- data
first_six$String <- substr(data$String, start = 1, stop = 6)

first_six <- first_six %>%
  group_by(String) %>%
  summarise(Sum = sum(Frequency))
