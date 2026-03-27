from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Union, overload, Iterable, Optional
import numpy as np 
from Crypto.Util.number import bytes_to_long, long_to_bytes 

# một số exception

class MatrixDimensionError(ValueError): 
    "Kích thước không hợp lệ"
    

class InvalidMatrixError(ValueError): 
    "Dữ liệu đầu vào không hợp lệ"
    

class InvalidParameterError(ValueError): 
    "Tham so khong hop le"
    
def mod_q(x: int, q: int) -> int:
    if q <= 0:
        raise InvalidParameterError
    return x % q 

def center_mod_q(x: int , q: int) -> int:
    if q <= 0: 
        raise InvalidParameterError
    x = x % q 
    if x >= q//2:
        return x - q
    else: 
        return x 

@dataclass
class Vector: 
    data: List[int] = field(default_factory = list)
    q : int = 0 

    def __post_init__(self): 
        self.data = [int(x) for x in self.data]
        if self.q > 0: 
            self.data = [mod_q(x, self.q) for x in self.data]
    
    @classmethod
    def from_iterable(cls, values: Iterable[int]) -> Vector:
        return cls(list(values))
    @classmethod 
    def zero(cls, size: int) -> Vector:
        return cls([0] * size)
    def __len__(self) -> int:
        return len(self.data) 
    def __iter__(self):
        return iter(self.data)
    def __getitem__(self, index: int) -> int: 
        return self.data[index]
    def __setitem__(self, index: int, value: int) -> None: 
        self.data[index] = value 
    def check_same_size(self, other: Vector) -> None: 
        if len(self) != len(other): 
            raise MatrixDimensionError
    def check_same_q(self, other: Vector) -> None: 
        if self.q != other.q: 
            raise InvalidParameterError
    # tinh toan 
    def __add__(self, other: Vector) -> Vector: 
        self.check_same_size(other) 
        self.check_same_q(other) 
        return Vector([a+b for a,b in zip(self.data, other.data)], q=self.q)
    def __sub__(self, other: Vector) -> Vector: 
        self.check_same_size(other) 
        self.check_same_q(other)
        return Vector([a-b for a,b in zip(self.data, other.data)], q=self.q)

    def __neg__(self) -> Vector:
        return Vector([-x for x in self.data], q=self.q)

    @overload
    def __mul__(self, other: int) -> Vector: ...
    @overload 
    def __mul__(self, other: Vector) -> int: ...
    def __mul__(self, other: Union[int, Vector]) -> Union[int, Vector]: 
        if isinstance(other, int): 
            return Vector([x*other for x in self.data], q=self.q)
        if isinstance(other, Vector): 
            self.check_same_size(other)
            self.check_same_q(other)
            return mod_q(sum(a*b for a,b in zip(self.data,other.data)), self.q)
        raise TypeError
    def __repr__(self) -> str: 
        return f"Vector({self.data})"
    def __rmul__(self, other: int) -> Vector: 
        if not isinstance(other, int): 
            raise TypeError
        
        return Vector([other *x for x in self.data], q=self.q) 
    def to_list(self)-> list[int]:
        return self.data.copy() 
    
# Matrix thi ta nen luu thong tin gi? Luu thong tin ve hang va cot cua no co can thiet hay k 
# Hay chi can luu Matrix duoi dang cac hang la duoc? 
@dataclass
class Matrix:
    rows : List[Union[List[int], Vector]] = field(default_factory=list)
    q : int = 0 

    def __post_init__(self) -> None: 
        if not self.rows: 
            raise InvalidMatrixError
        
        self.rows = [Vector(r, q=self.q) if not isinstance(r,Vector) else r for r in self.rows]
        row_len = [len(x) for x in self.rows] 
        if row_len[0] == 0:
            raise InvalidMatrixError
        if any(length != row_len[0] for length in row_len):
            raise InvalidMatrixError
        
    @classmethod
    def zero(cls,nrows:int, ncols:int) -> Matrix: 
        res = [[0 for _ in range(ncols)] for _ in range(nrows)]
        return cls(res)
    @property
    def nrows(self)-> int: 
        return len(self.rows)
    @property 
    def ncols(self) -> int:
        return len(self.rows[0])
    def __repr__(self): 
        return "Matrix([\n  " + ",\n  ".join([str(v.data) for v in self.rows]) + "\n])"
    def __getitem__(self, index: int) -> Vector: 
        return self.rows[index]
    def check_same_size(self, other: Matrix) -> None: 
        if self.nrows != other.nrows or self.ncols != other.ncols: 
            raise MatrixDimensionError
    def check_same_q(self, other: Matrix) -> None:
        if self.q != other.q: 
            raise InvalidParameterError 
    
    
    def __add__(self, other: Matrix) -> Matrix: 
        self.check_same_size(other) 
        self.check_same_q(other)
        return Matrix([r1+r2 for r1,r2 in zip(self.rows, other.rows)],q=self.q)
    def __sub__(self, other: Matrix) -> Matrix: 
        self.check_same_size(other)
        self.check_same_q(other)
        return Matrix([r1 - r2 for r1,r2 in zip(self.rows, other.rows)],q=self.q)
    def __neg__(self) -> Matrix: 
        return Matrix([-r for r in self.rows], q=self.q)

    def __mul__(self, other: int) -> Matrix:  
        if isinstance(other, int): 
            return Matrix([[x*other for x in r] for r in self.rows], q=self.q)
        raise TypeError
    def __matmul__(self, other: Union[Vector, Matrix]) -> Union[Vector, Matrix]:
        if isinstance(other, Vector): 
            if self.q != other.q:
                raise InvalidParameterError
            if self.ncols != len(other):
                raise MatrixDimensionError
            res = [sum(self.rows[i][k] * other[k] for k in range(self.ncols)) for i in range(self.nrows)]
            return Vector(res, q=self.q)
        if isinstance(other, Matrix): 
            if self.ncols != other.nrows: 
                raise MatrixDimensionError 
            if self.q != other.q: 
                raise InvalidParameterError
            res = [[0 for _ in range(other.ncols)] for _ in range(self.nrows)]
            for i in range(self.nrows): 
                for j in range(other.ncols): 
                    res[i][j] = sum(self.rows[i][k] * other[k][j] for k in range(self.ncols))
            return Matrix(res, q=self.q)
        
    def row(self, index: int): 
        return self.rows[index]
    def col(self, index: int):
        return Vector([r[index] for r in self.rows], q = self.q)



INT32_MIN = np.iinfo(np.int32).min
INT32_MAX = np.iinfo(np.int32).max

def uniform_sample_int32(size: int , q: int) -> Vector: 
    return Vector(np.random.randint(INT32_MIN, INT32_MAX+1, size, dtype = np.int32).tolist(), q=q)

def gaussian_sample_int32(std: float, size: Optional[int]) -> Vector:
    if size is None:
        size = 1
    samples = INT32_MAX * np.random.normal(0.0, std, size)

    samples = np.clip(
        samples,
        INT32_MIN,
        INT32_MAX
    )

    return Vector(np.round(samples).astype(np.int32).tolist(), q=0)
def gaussian_sample_int32_scalar(std: float) -> int:
    return gaussian_sample_int32(std, size=None)[0]

# LWE scheme 


@dataclass(frozen=True)
class LWEParams:
    # SIZE OF ENCRYPTION KEY 
    n: int 
    q : int 
    noise_bound: int 



@dataclass
class LWEPlaintext:
    message: int 


@dataclass
class LWECiphertext:
    params: LWEParams
    a: Vector 
    b: int 
    # Z^n_q x Z_qq


@dataclass
class LWEKey: 
    params: LWEParams
    key: Vector 


def generate_lwe_key(params: LWEParams) -> LWEKey: 
    key = Vector(np.random.randint(low = 0 , high = 2, size = params.n, dtype = np.int32).tolist(), q=params.q)
    return LWEKey(params, key)

def lwe_encrypt(plaintext: LWEPlaintext, key: LWEKey) -> LWECiphertext: 
    params = key.params 
    e = gaussian_sample_int32_scalar(std = params.noise_bound) 

    a = uniform_sample_int32(params.n, q=params.q)
    b = (a*key.key + e + plaintext.message) % params.q 
    return LWECiphertext(params, a,b)


def lwe_decrypt(ciphertext: LWECiphertext, key: LWEKey) -> LWEPlaintext:
    params = key.params
    a,b = ciphertext.a, ciphertext.b 
    pt = (b - a*key.key) % params.q 
    return LWEPlaintext(pt)


def encode(m : int) -> int: 
    return m*(2**29)
def decode(m: int) -> int: 
    return int(round(m/(2**29))) % 8

def LWE_Encode(i: int) -> LWEPlaintext: 
    return LWEPlaintext(encode(i))

def LWE_Decode(pt: LWEPlaintext) -> int:
    return decode(pt.message)




param = LWEParams(n=1024, q = 2**32, noise_bound = 2**(-24)) 


key = generate_lwe_key(param)
plaintext = LWE_Encode(4) 
ciphertext = lwe_encrypt(plaintext, key)
dec_ciphertext = lwe_decrypt(ciphertext, key)
decoded = LWE_Decode(dec_ciphertext)
print(decoded) 
assert decoded == 4, "Decryption failed"



# test a msg
from Crypto.Util.number import bytes_to_long, long_to_bytes 

msg = b'CCTF{s0lV1n9_4_Syst3m_0f_L1n3Ar_3qUaTi0n5_0vEr_7H3_F!3lD_F(2)!}'

# chuyen moi bytes -> octet

def bytes_to_octet(msg: bytes) -> list[int]: 
    digits = []
    for b in msg: 
        oct = format(b, '03o')
        digits.extend(int(d) for d in oct) 
    return digits 

def octet_to_bytes(octets: list[int]) -> bytes: 
    if len(octets) % 3 != 0:
        raise ValueError("Invalid octet length")
    out = bytearray() 
    for i in range(0, len(octets), 3): 
        d0,d1,d2 = octets[i], octets[i+1],octets[i+2]
        if not (0<=d0 <8 and 0<=d1<8 and 0<=d2<8): 
            raise ValueError("Invalid digits") 
        value = d0*64+d1*8+d2
        out.append(value) 
    return bytes(out) 




def lwe_encrypt_message(msg: bytes, key: LWEKey) -> list[LWECiphertext]: 
    m = bytes_to_octet(msg) 
    for i in range(len(m)): 
        m[i] = LWE_Encode(m[i]).message
    ciphertexts = [lwe_encrypt(LWEPlaintext(m[i]), key) for i in range(len(m))]
    return ciphertexts

def lwe_decrypt_message(ciphertexts: list[LWECiphertext], key: LWEKey) -> bytes: 
    plaintexts = [lwe_decrypt(ct, key) for ct in ciphertexts]
    octets = [LWE_Decode(pt) for pt in plaintexts]
    return octet_to_bytes(octets)
ct = lwe_encrypt_message(msg, key)
decrypted_msg = lwe_decrypt_message(ct, key)
print(decrypted_msg)
# if __name__ == "__main__": 
#     pass 


