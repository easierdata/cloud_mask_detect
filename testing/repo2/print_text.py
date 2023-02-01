def print_text():
    with open("inputs/hello.txt", "r") as f:
        text = f.read()
        print(text)


if __name__ == "__main__":
   print_text() 