import os
input_dir = "inputs/"
# Get the list of all files in the directory
dir_list = os.listdir(input_dir)

# Iterate through the list and find the first .txt file
for file in dir_list:
    if file.endswith(".txt"):
        file_path = os.path.join(input_dir, file)
        with open(file_path, "r") as f:
            contents = f.read()
            with open("outputs/container_output.txt", "w") as out:
                out.write(contents+" - This is a test")
                break

# Print the contents of the output file
with open("outputs/container_output.txt", "r") as f:
    print(f.read())