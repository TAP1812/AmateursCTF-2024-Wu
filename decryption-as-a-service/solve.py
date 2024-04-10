
from pwn import *
from Cryptodome.Util.number import isPrime, long_to_bytes

HOST = 'chal.amt.rs'
PORT = 1417

primes = []
for i in range(5000):
    if isPrime(i):
        primes.append(i)

while True:
    x = []
    r = remote(HOST, PORT)
    c = int(r.recvline().decode().strip().split("encrypted_flag = ")[1])

    for i in range(4):
        r.sendlineafter(b'message? ', str(pow(2,1024 + i)).encode())
        x.append(int(r.recvline().decode().strip(),16))

    N = math.gcd(x[0] * x[2] - x[1] * x[1], x[1] * x[3] - x[2] * x[2])
    if c % 2 ==0:
        print(f"c = {c}")
        for prime in primes:
            while N % prime ==0:
                N = N // prime

        try:
            first_part = (x[1] * pow(x[0], -1, N)) % N
            r.sendlineafter(b'message? ', str(c // 2).encode())
            second_part = int(r.recvline().decode().strip(),16)
            plain_text = (first_part * second_part) % N
            print(long_to_bytes(plain_text).decode())
            r.close()
            break
        except:
            r.close()
    else:
        r.close()



