##### Code to change xlsx > tsv

import pandas as pd

# Define the input and output file paths
input_file = r"paste path to the xlsx file directory here.xlsx"
output_file = r"paste path to the newly created tsv file directory here.tsv"


# Read the Excel file
df = pd.read_excel(input_file)  # By default, reads the first sheet

# Save the DataFrame to a .tsv file
df.to_csv(output_file, sep="\t", index=False)

print(f"Conversion complete! Saved as {output_file}.")
