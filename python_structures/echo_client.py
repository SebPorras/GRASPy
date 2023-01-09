############################################
# Date: 2/12/22
# Author: Sebastian Porras
# Aims: Create a simple echo practice server
# and client to practice socket programming
# # Based off tutorial https://realpython.com/python-sockets/#reference
############################################

import socket
import message


test = message.constructJsonMessage(
    aln="./small_test_data/test_aln.aln",
    nwk="./small_test_data/test_nwk.nwk",
    type="Joint"
)

HOST = "127.0.0.1"
PORT = 4072

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))

    # convert to bytes and send
    s.sendall(test.encode())

    data = s.recv(1024)

print(f"Response: {data!r}")
