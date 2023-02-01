import os
input_dir = "inputs/LC08_L1TP_152028_20160209_20200907_02_T1/"
# Get the list of all files in the directory
dir_list = os.listdir(input_dir)

# Iterate through the list and find the first .txt file
for file in dir_list:
    if file.endswith(".txt"):
        file_path = os.path.join(input_dir, file)
        with open(file_path, "r") as f:
            contents = f.read()
            with open("outputs/container_output.txt", "w") as out:
                out.write(contents)
                break

# Print the contents of the output file
with open("outputs/container_output.txt", "r") as f:
    # print first 5 lines
    for i in range(5):
        print(f.readline())