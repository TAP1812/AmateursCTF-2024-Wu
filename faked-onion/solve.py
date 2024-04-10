from pwn import remote
from string import printable
from tqdm import tqdm

def strxor(a: bytes, b: bytes):
    return bytes([x ^ y for x, y in zip(a, b)])

HOST = 'chal.amt.rs'
PORT = 1414

r = remote(HOST, PORT)
r.recvline()

r.sendlineafter(b'3. Exit\n> ', b'2')
enc = bytes.fromhex(r.recvline().decode().strip())

hash = dict()
for i in tqdm(range(100)):
    payload = printable[i] * 16
    r.sendlineafter(b'3. Exit\n> ', b'1')
    r.sendlineafter(b'Enter your message in hex: ', payload.encode().hex().encode())
    hashed = strxor(bytes.fromhex(r.recvline().decode().strip())[1:], payload[:15].encode())
    hash[printable[i]] = hashed

flag = b''
for i in range(0, len(enc), 16):
    block = enc[i : i + 16]
    flag += strxor(block[1:], hash[chr(block[0])]) + chr(block[0]).encode()
    print(flag)
