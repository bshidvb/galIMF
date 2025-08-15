import re

# File paths
input_file = "/Users/adriana_work/Desktop/galIMF/yield_tables__2024/TABLE1-VW93ML_total_yields_scaled_with_X_i_0.txt"
output_file = "/Users/adriana_work/Desktop/galIMF/yield_tables__2024/TABLE1-VW93ML_total_yields_scaled_with_X_i_0_total_yield.txt"

# Regular expression to match rows starting with '&'
row_pattern = re.compile(r"^&")

# Read the input file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Process the file
updated_lines = []
skipped_rows = []
for line in lines:
    if row_pattern.match(line):
        # Split the line into columns
        columns = line.split()

        # Extract relevant values
        try:
            p_i = float(columns[1])  # Net Yield
            M_ej = float(columns[2])  # MassExp
            X_i_0 = float(columns[-1])  # X_i,0 (last column)

            # Calculate Total yield
            total_yield = p_i + M_ej * X_i_0

            # Update the Total yield column (8th column, index 7)
            columns[7] = f"{total_yield:.6e}"
        except (ValueError, IndexError) as e:
            # Log skipped rows for debugging
            skipped_rows.append((line.strip(), str(e)))
            continue

        # Reconstruct the line
        line = " ".join(columns) + "\n"

    # Append the processed or unprocessed line
    updated_lines.append(line)

# Write the updated lines back to the file
with open(output_file, 'w') as file:
    file.writelines(updated_lines)

# Log skipped rows
if skipped_rows:
    print("Skipped rows:")
    for row, error in skipped_rows:
        print(f"Row: {row}, Error: {error}")
else:
    print("All rows processed successfully.")

print("Total yield column updated successfully.")
