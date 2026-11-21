# extract per cell type lifespan (days) from Sender and Milo 2021 Summary.xlsx; rodent value when human is missing
import openpyxl, sys
wb = openpyxl.load_workbook(sys.argv[1], data_only=True)
with open(sys.argv[2], "w") as f:
    f.write("cell_type\tlifespan_days\tlifespan_rodent_days\n")
    for ws in wb.worksheets[1:]:
        d = {r[0]: r[1] for r in ws.iter_rows(min_row=2, values_only=True) if r[0]}
        f.write(f"{ws.title}\t{d.get('lifespan')}\t{d.get('lifespan in rodents')}\n")
