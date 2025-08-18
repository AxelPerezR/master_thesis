import pandas as pd
import matplotlib.pyplot as plt

mapping_file = "protoboard2crystal mapping.csv" # Path to the mapping file # Provided by Olivier
rms_file = "RMS_values_fpga2.csv" # Path to the RMS per channel file
pb = "06" # Protoboard serial number # fpga0=01, fpga1=008, fpga2=06

connectors = {
    "A": "JA pin",
    "B": "JB pin",
    "C": "JC pin",
    "D": "JD pin"
}

# --- Load and clean mapping file ---
mapping_df = pd.read_csv(mapping_file, delimiter=';')
mapping_df.columns = [col.strip() for col in mapping_df.columns]

# --- Load and prepare RMS file ---
rms_df = pd.read_csv(rms_file)

# Create subplot grid
fig, axs = plt.subplots(2, 4, figsize=(16, 8), sharey=True)
colors = ['red', 'blue', 'green']
labels = ['No cable connected', 'Cable connected', 'Cable + Voltage']

for i, asic in enumerate([0, 1]):  # Rows
    for j, (letter, connector) in enumerate(connectors.items()):  # Columns
        # --- Prepare mapping data ---
        map_df = mapping_df[['PB_INA', connector]].dropna()
        map_df['PB_INA'] = map_df['PB_INA'].astype(int)
        map_df[connector] = map_df[connector].astype(int)

        # Remove duplicate PB_INA according to ASIC
        _keep = 'first' if asic == 0 else 'last'
        rms_filtered = rms_df.drop_duplicates(subset='PB_INA', keep=_keep)

        # Merge on PB_INA
        merged = pd.merge(map_df, rms_filtered[['PB_INA', 'rms_no', 'rms_cab', 'rms_cab_hv']], on='PB_INA', how='left')

        # Sort connector pins: odd first, then even
        odd = merged[merged[connector] % 2 == 1].sort_values(connector)
        even = merged[merged[connector] % 2 == 0].sort_values(connector)
        sorted_df = pd.concat([odd, even])

        # --- Plotting ---
        ax = axs[i, j]
        x_vals = sorted_df[connector].astype(str)

        ax.scatter(x_vals, sorted_df['rms_no'], label=labels[0], color=colors[0])
        ax.scatter(x_vals, sorted_df['rms_cab'], label=labels[1], color=colors[1])
        ax.scatter(x_vals, sorted_df['rms_cab_hv'], label=labels[2], color=colors[2])

        ax.set_title(f'ASIC {asic} - Cable {letter}')
        ax.set_xlabel(f'J{letter} pin (odd→even)')
        if j == 0:
            ax.set_ylabel('RMS (ADC)')
        ax.set_ylim(0, 30)
        ax.grid(True)
        if i == 0 and j == 0:
            ax.legend(loc=1)

plt.suptitle(f'RMS Noise Measurements by ASIC and Connector (Protoboard {pb})', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
