############################################
# Date: 2/12
# Author: Sebastian Porras
# Aims: Create a simple echo practice server
# and client to practice socket programming
# Based off tutorial https://realpython.com/python-sockets/#reference
############################################

import socket

HOST = "127.0.0.1"
PORT = 4072

# af_inet is internet address family for IPv4
# sock stream is the TCP (protocol) for network
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:

    dat = ""
    # PORT is the number that accepts connections from the client

    s.bind((HOST, PORT))
    s.listen()

    # new socket object is contained in conn
    conn, addr = s.accept()

    with conn:
        print(f"Connected by {addr}")

        i = 0
        while True:

            # recieves data packages
            data = conn.recv(1024)
            dat += data.decode()

            if not data:
                break
            i += 1
            conn.sendall(f'Recieved package {i}, '.encode())

    print(f"Data recieved {dat}")
