###############################################################################
# Date: 20/1/23
# Author: Sebastian Porras
# Aims: A simple client socket that allows the user to send and recieve
# JSON strings.
###############################################################################

import socket
import time
import sys


def receive_message(socket, timeout=2) -> str:
    """Recieves chunks of data from the server 
    and decodes this from bytes back into a 
    string. 
    """

    # allows the socket to recieve multiple chunks
    socket.setblocking(False)

    # chunks will be appended here
    message = []
    chunk = ''

    start = time.time()

    while True:

        # if you have data, break after the timeout
        if message and time.time() - start > timeout:

            break

        # check for when no data is coming in
        elif time.time() - start > timeout*2:

            break

        try:
            # read in the data and record it
            chunk = socket.recv(8192)
            if chunk:
                toStr = chunk.decode()
                message.append(toStr)

                # reset the timeout
                start = time.time()

            else:
                time.sleep(0.1)

        # if no data, continue and re-test conditions
        except:
            pass

    return ''.join(message)


def sendRequest(message: str) -> str:
    """User enters their message which is 
    converted into bytes before being sent to 
    the server. Also recieves the response and 
    returns this to the user. 
    """

    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    except socket.error:

        print("Socket failed to be created")
        sys.exit()

    #print("Socket created...\n")

    HOST = '10.139.1.21'
    PORT = 4072

    print("Connecting to server...\n")

    s.connect((HOST, PORT))

    print(f"Socket connected to {HOST} on IP {PORT}\n")

    try:
        s.sendall(message.encode())

    except socket.error:

        print("Send failed\n")
        sys.exit()

    response = receive_message(s)

    #print("Closing socket...\n")

    s.close()

    return response
