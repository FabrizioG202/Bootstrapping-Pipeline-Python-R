def red(text : object) -> str:
    return "\033[31m" + str(text) + "\033[0m"

#the blue color
def blue(text : object) -> str:
    return "\033[34m" + str(text) + "\033[0m"

#the green color
def green(text : object) -> str:
    return "\033[32m" + str(text) + "\033[0m"

#the yellow color
def yellow(text : object) -> str:
    return "\033[33m" + str(text) + "\033[0m"

#the cyan color
def cyan(text : object) -> str:
    return "\033[36m" + str(text) + "\033[0m"

# the magenta color
def magenta(text : object) -> str:
    return "\033[35m" + str(text) + "\033[0m"

# green background
def green_background(text : object) -> str:
    return "\033[42m" + str(text) + "\033[0m"

# red background
def red_background(text : object) -> str:
    return "\033[41m" + str(text) + "\033[0m"

# bold
def bold(text : object) -> str:
    return "\033[1m" + str(text) + "\033[0m"

# underline
def underline(text : object) -> str:
    return "\033[4m" + str(text) + "\033[0m"

def italic(text : object) -> str:
    return "\033[3m" + str(text) + "\033[0m"

# white background
def white_background(text : object) -> str:
    return "\033[47m" + str(text) + "\033[0m"

