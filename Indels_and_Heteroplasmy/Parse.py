# Author: Miguel Campos
# Date: January 2025

import pandas as pd

# Load the indels and heteroplasmic positions files
indels_path = 'Indels.csv'
heteroplasmic_positions_path = 'Heteroplasmy.csv'

# Read the indels file
indels = pd.read_csv(indels_path, delimiter='\t')

# Read the heteroplasmic positions file
heteroplasmic_positions = pd.read_csv(heteroplasmic_positions_path, delimiter='\t')

# Create a list to store the results
results = []

# Iterate through each heteroplasmic position
for idx, heteroplasmy in heteroplasmic_positions.iterrows():
    position = heteroplasmy['Position']
    # Check if the position is within any indel region
    for _, indel in indels.iterrows():
        if indel['Start'] <= position <= indel['End']:
            results.append({
                'Position': position,
                'B.stacei': heteroplasmy['B.stacei'],
                'B.hybridumS': heteroplasmy['B.hybridumS'],
                'B.hybridumD': heteroplasmy['B.hybridumD'],
                'B.distachyon': heteroplasmy['B.distachyon'],
                'Indel_Start': indel['Start'],
                'Indel_End': indel['End'],
                'Indel_Size': indel['Size'],
                'Indel_Type': indel['Type']
            })

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)

# Save the results to a new CSV file
results_df.to_csv('heteroplasmy_in_indels.csv', index=False)

print("Results saved in heteroplasmy_in_indels.csv")
