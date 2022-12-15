import tkinter as tk
import tkinter.filedialog
from excel_parser import parse_excel_input

# Open a file browser and select the input file
def select_input_file():
    root = tk.Tk()
    root.withdraw()
    input_file = tkinter.filedialog.askopenfilename(
        title="Select the input Excel file",
        filetypes=[("Excel files", ".xlsx .xlsm")],
    )
    return input_file


def main():
    input_file = select_input_file()
    print(input_file)
    parse_excel_input(input_file)


if __name__ == "__main__":
    main()
