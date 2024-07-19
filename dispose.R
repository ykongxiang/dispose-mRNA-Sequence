#mrna_file 要包含Sequence,Strand,FS_start,FS_end列 且列含义相对应
if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
library(Biostrings)

extract_sequences <- function(mrna_file) {
  mrna_data <- read_excel(mrna_file)
  mrna_data <- na.omit(mrna_data)  
  mrna_data$FS_start <- as.numeric(mrna_data$FS_start)
  mrna_data$FS_end <- as.numeric(mrna_data$FS_end)
  
  extract_and_reverse_substring <- function(text, start, end, strand) {
    if (strand == '-') {
      text <- as.character(reverseComplement(DNAString(text)))
    }
    # 提取子串
    substring <- substr(text, start - 2, end + 1)
    if (start > end){
      substring <- substr(text,end - 2, start + 1)
    }
    return(substring)
  }
  
  seqs <- sapply(1:nrow(mrna_data), function(i) {
    extract_and_reverse_substring(
      text = mrna_data$Sequence[i],
      start = mrna_data$FS_start[i],
      end = mrna_data$FS_end[i],
      strand = mrna_data$Strand[i]
    )
  })
  return(seqs)
}
#####调用示例：
set_data_frame <- function(x){
  string_counts <- table(x)
  total_count <- sum(string_counts)
  string_frequencies <- string_counts / total_count
  data <- data.frame(
    String = names(string_frequencies),
    Frequency = as.numeric(string_frequencies)
  )
  data <- data[order(-data$Frequency), ]
  return(data)
}
a <- set_data_frame(extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx'))

first_six <-function(data){
  first_six$String <- substr(data$String, start = 1, stop = 6)

  first_six <- first_six %>%
    group_by(String) %>%
    summarise(Sum = sum(Frequency))
  return(first_six)
}
