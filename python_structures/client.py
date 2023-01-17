############################################
# Date: 2/12/22
# Author: Sebastian Porras
# Aims: Create a simple echo practice server
# and client to practice socket programming
# # Based off tutorial https://realpython.com/python-sockets/#reference
############################################

import socket


def sendRequest(message: str) -> str:
    """
    Constructs a socket and sends a message to the 
    bnkit sever

    Paramters:
        message(str): Any generic message

    Returns:
        Response from bnkit server
    """

    #HOST = "127.0.0.1"
    HOST = '10.139.1.21'

    PORT = 4072

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM)\
            as s:

        print("Establishing connection...\n")
        s.connect((HOST, PORT))

        s.sendall(message.encode())
        print("Request sent...\n")

        response = s.recv(2048)

    print("Response ")
    print(f"{response!r}")
    return f"Response: {response!r}"
