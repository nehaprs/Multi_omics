# Load necessary libraries
library(readxl)
library(openxlsx)
library(readxl)
library(dplyr)
library(writexl)
rna<- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")

rna1 = rna[,c(1,14,17)]
rna1$FC = 2^(rna1$logFC)
colnames(rna1) = c("Protein", "Transcriptome log FC", "Transcriptome FDR", "Transcriptome FC")

PTT_tableS2 <- read_excel("PTT_tableS2.xlsx")
PTT = inner_join(PTT_tableS2, rna1, by = "Protein")
write_xlsx(PTT,"PTT_tableS2v2.xlsx")

TT_tableS1 <- read_excel("TT_tableS1.xlsx")
TT = inner_join(TT_tableS1, rna1, by = "Protein")




# Define a function to process each xlsx file
process_excel_file <- function(input_file) {
  
  # Read data from Sheet1 and Sheet2 of the input file
  sheet1_data <- read_excel(input_file, sheet = 1)
  sheet2_data <- read_excel(input_file, sheet = 2)
  
  # Perform operations on Sheet1 and Sheet2
  
    result_sheet1 = inner_join(sheet1_data, rna1, by = "Protein")
    result_sheet2 = inner_join(sheet2_data, rna1, by = "Protein")
  
  # Convert the result back to a data frame
  result_sheet1 <- as.data.frame(result_sheet1) 
  result_sheet2 <- as.data.frame(result_sheet2)  
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add the first sheet and write the result of Sheet1 operation
  addWorksheet(wb, "Sheet1")
  writeData(wb, sheet = 1, result_sheet1)
  
  # Add the second sheet and write the result of Sheet2 operation
  addWorksheet(wb, "Sheet2")
  writeData(wb, sheet = 2, result_sheet2)
  output_file = paste0(input_file,"v2")
  # Save the new workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Print a message indicating that processing is done
  print(paste("Results written to", output_file))
}

process_excel_file("TT_tableS1.xlsx")


###########################################################

#annotation
pttDEP <- read_excel("~/BINF/BINF_old/multi-omics/transcriptome-cp redo/analysis after integration/pttDEP.xlsx", 
                     +     sheet = "Sheet2")
pttDEP2 = inner_join(pttDEP, rna1, by = "Protein")
pttall = as.data.frame(pttDEP2)



# Define a function to process each xlsx file
process_annotation <- function(input_file) {
  
  # Read data from Sheet1 and Sheet2 of the input file
  for(i in 1:7){
    print(i)
    sheetdata = as.data.frame(read_excel(input_file, sheet = i))
    print(head(sheetdata))
    
    resultsheet = inner_join(sheetdata, pttall, by = "Protein")
    paste0("result_sheet",i) = resultsheet
    
  }
  
  
  # Create a new workbook
  wb <- createWorkbook()
  # Add the sheets
  for(i in 1:7){
    sheeti = paste0("Sheet",i)
    addWorksheet(wb, sheeti)
    writeData(wb, sheet = i, paste0("result_sheet",i) )
  }
  
  
    output_file = paste0(input_file,"v2")
  # Save the new workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Print a message indicating that processing is done
  print(paste("Results written to", output_file))
}


# Load necessary libraries
library(readxl)
library(openxlsx)
library(dplyr)
############################################
# Define a function to process each xlsx file
process_annotation <- function(input_file, pttall) {
  
  # Initialize a list to store the result for each sheet
  result_list <- list()
  
  # Read data from Sheet1 to Sheet7 of the input file
  for (i in 1:7) {
    print(paste("Processing sheet", i))
    
    # Read the current sheet as a data frame
    sheetdata <- as.data.frame(read_excel(input_file, sheet = i))
    
    
    # Perform the join operation with 'pttall'
    resultsheet <- inner_join(sheetdata, pttall, by = "Protein")
    
    # Store the result in the list with the name result_sheet_i
    result_list[[i]] <- resultsheet
  }
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add the sheets to the workbook and write the result data from the list
  for (i in 1:7) {
    sheeti <- paste0("Sheet", i)
    addWorksheet(wb, sheeti)
    
    # Write the corresponding result from the list to the new sheet
    writeData(wb, sheet = i, result_list[[i]])
  }
  
  # Define the output file name
  output_file <- paste0(input_file, "_v2.xlsx")
  
  # Save the new workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Print a message indicating that processing is done
  print(paste("Results written to", output_file))
}

###################
# process_annotation("your_file.xlsx", pttall)


#apply
process_annotation("annotation_TableS5.xlsx", pttall)
