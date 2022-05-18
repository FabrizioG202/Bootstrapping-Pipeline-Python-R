def red(text):
    return "\033[91m" + text + "\033[0m"

def green(text):
    return "\033[92m" + text + "\033[0m"

def blue(text):
    return "\033[94m" + text + "\033[0m"

def log_info(text):
    print(blue(text))

def yellow(text):
    return "\033[93m" + text + "\033[0m"

def bold(text):
    return "\033[1m" + text + "\033[0m"