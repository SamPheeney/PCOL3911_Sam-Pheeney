# this line installs the necessary libraries
# it gets the list of libraries from the file requirements.txt
python3 -m pip install -r requirements.txt
#this line runs the python file
python3 generate_fp_from_smiles.py
data = pd.read_csv("gsar_a_1049665_sm4518.csv")
