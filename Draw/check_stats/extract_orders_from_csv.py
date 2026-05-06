import json
import pandas as pd
from prettytable import PrettyTable

# Cross-sections in pb for the different jet multiplicities
XS0=1788
XS1=339.7
XS2=125.1
XS = [XS0, XS1, XS2]
XS_tot = sum(XS)

# Load the CSV file into a DataFrame
csv_file_path = 'Draw/check_stats/DYto2Tau.csv'
df = pd.read_csv(csv_file_path)

columns_of_interest = [
    'Dataset', 
    'Root request status',
    'Chained request',
    'NanoAOD done events'
]

# Filter the DataFrame to include only the columns of interest
df = df[columns_of_interest]
mask = df['Root request status'] == 'done'
df = df[mask]

# Separate the DataFrame into different years
df_2022 = df[df['Chained request'].str.contains('Run3Summer22NanoAODv12')]
df_2022EE = df[df['Chained request'].str.contains('Run3Summer22EENanoAODv12')]
df_2023 = df[df['Chained request'].str.contains('Run3Summer23NanoAODv12')]
df_2023BPix = df[df['Chained request'].str.contains('Run3Summer23BPixNanoAODv12')]
df_2024 = df[df['Chained request'].str.contains('RunIII2024Summer24NanoAODv15')]

# Dictionary to be populated
dictionary = {
    'Run3_2022': {
        'DY_to2Tau_0J': int(df_2022['NanoAOD done events'].iloc[3]),
        'DY_to2Tau_1J': int(df_2022['NanoAOD done events'].iloc[4]),
        'DY_to2Tau_2J': int(df_2022['NanoAOD done events'].iloc[5]),
        'DY_to2Tau_0J_filtered': int(df_2022['NanoAOD done events'].iloc[0]),
        'DY_to2Tau_1J_filtered': int(df_2022['NanoAOD done events'].iloc[1]),
        'DY_to2Tau_2J_filtered': int(df_2022['NanoAOD done events'].iloc[2]), 
    },
    'Run3_2022EE': {
        'DY_to2Tau_0J': int(df_2022EE['NanoAOD done events'].iloc[3]),
        'DY_to2Tau_1J': int(df_2022EE['NanoAOD done events'].iloc[4]),
        'DY_to2Tau_2J': int(df_2022EE['NanoAOD done events'].iloc[5]),
        'DY_to2Tau_0J_filtered': int(df_2022EE['NanoAOD done events'].iloc[0]),
        'DY_to2Tau_1J_filtered': int(df_2022EE['NanoAOD done events'].iloc[1]),
        'DY_to2Tau_2J_filtered': int(df_2022EE['NanoAOD done events'].iloc[2]), 
    },
    'Run3_2023': {
        'DY_to2Tau_0J': int(df_2023['NanoAOD done events'].iloc[3]),
        'DY_to2Tau_1J': int(df_2023['NanoAOD done events'].iloc[4]),
        'DY_to2Tau_2J': int(df_2023['NanoAOD done events'].iloc[5]),
        'DY_to2Tau_0J_filtered': int(df_2023['NanoAOD done events'].iloc[0]),
        'DY_to2Tau_1J_filtered': int(df_2023['NanoAOD done events'].iloc[1]),
        'DY_to2Tau_2J_filtered': int(df_2023['NanoAOD done events'].iloc[2]), 
    },
    'Run3_2023BPix': {
        'DY_to2Tau_0J': int(df_2023BPix['NanoAOD done events'].iloc[3]),
        'DY_to2Tau_1J': int(df_2023BPix['NanoAOD done events'].iloc[4]),
        'DY_to2Tau_2J': int(df_2023BPix['NanoAOD done events'].iloc[5]),
        'DY_to2Tau_0J_filtered': int(df_2023BPix['NanoAOD done events'].iloc[0]),
        'DY_to2Tau_1J_filtered': int(df_2023BPix['NanoAOD done events'].iloc[1]),
        'DY_to2Tau_2J_filtered': int(df_2023BPix['NanoAOD done events'].iloc[2]),
    },
    'Early_Run3': {},
    'Run3_2024': {
        'DY_to2Tau_0J': int(df_2024['NanoAOD done events'].iloc[0]),
        'DY_to2Tau_1J': int(df_2024['NanoAOD done events'].iloc[1]),
        'DY_to2Tau_2J': int(df_2024['NanoAOD done events'].iloc[2]),
        'DY_to2Tau': int(df_2024['NanoAOD done events'].iloc[4]),
    },
}

early_run3_eras = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
for sample in dictionary['Run3_2022'].keys():
    events = sum(dictionary[era][sample] for era in early_run3_eras)
    dictionary['Early_Run3'][sample] = events

table = PrettyTable()
table.title = "Number of DYto2Tau events ordered"
table.field_names = ["n_jets", "Early_Run3", "Run3_2024"]
for i, NJ in enumerate(['0J', '1J', '2J']):
    table.add_row([
        NJ,
        f"{dictionary['Early_Run3'][f'DY_to2Tau_{NJ}'] + dictionary['Early_Run3'][f'DY_to2Tau_{NJ}_filtered']:.2e}",
        f"{dictionary['Run3_2024'][f'DY_to2Tau_{NJ}'] + (XS[i] / XS_tot) * dictionary['Run3_2024'][f'DY_to2Tau']:.2e}"
    ])
print(table)

# Save the dictionary to a JSON file
json_file_path = 'Draw/check_stats/DYto2Tau.json'
with open(json_file_path, 'w') as json_file:
    json.dump(dictionary, json_file, indent=4)