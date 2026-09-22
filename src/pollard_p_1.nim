# Integer factorization via Pollard's P-1 method
#
# Note: gcd is the one provided by the bigints package.

import bigints, algorithm, sequtils, trutils

const
  zero = 0.initBigInt
  one = 1.initBigInt
  two = 2.initBigInt
  three = 3.initBi
  five = 5.initBigInt
  six = 6.initBigInt

# Return all primes up to n (sieve of Eratosthenes)
func primesUpTo(n: int): seq[int] =
  if n < 2:
    return @[]
  var sieve = newSeq[bool](n + 1)
  for i in 0 .. n:
    sieve[i] = true
  sieve[0] = false
  sieve[1] = false
  var i = 2
  while i * i <= n:
    if sieve[i]:
      for j in countup(i * i, n, i):
        sieve[j] = false
    i += 1
  for i in 2 .. n:
    if sieve[i]:
      result.add(i)

# Simple primality test (for small numbers)
func isPrime(n: BigInt): bool =
  if n < two:
    return false
  if n == two or n == three:
    return true
  if n mod two == zero or n mod three == zero:
    return false
  var i = five
  while i * i <= n:
    if n mod i == zero or n mod (i + two) == zero:
      return false
    i = i + six
  return true

# Integer factorization via Pollard's P-1 method
#
# n: positive integer to factor
# B1: bound for stage 1 (small primes)
# B2: bound for stage 2 (slightly larger primes)
# a: base (usually 2)
#
# Return a factor of n (or 1 if none is found). If n is prime, return n itself.
func pollardP1(n: BigInt, B1: int = 1000, B2: int = 100000, a: BigInt = two): BigInt =
  # Check whether n is prime
  if n < two:
    return one
  if isPrime(n):
    return n  # return n itself (prime)

  # Check that a and n are coprime
  let g = gcd(a, n)
  if g != one:
    return g

  # ---- Stage 1: for each prime p <= B1, compute x = x^p mod n iteratively ----
  # Compute step by step instead of one huge exponentiation (more memory-efficient)
  var x = a
  for p in primesUpTo(B1):
    x = powmod(x, p.initBigInt, n)

  # Check whether stage 1 found a factor
  let d1 = gcd(x - one, n)
  if one < d1 and d1 < n:
    return d1
  if d1 == n:
    # failure
    return one

  # ---- Stage 2: process primes B1 < q <= B2 ----
  for q in primesUpTo(B2):
    if q <= B1:
      continue
    # x = x^q mod n
    x = powmod(x, q.initBigInt, n)
    # Check gcd at each step (costly, but allows early detection)
    # Here we simply check at every step for brevity
    let d = gcd(x - one, n)
    if one < d and d < n:
      return d
    if d == n:
      # failed in stage 2 as well
      break

  return one

# Factorize by repeatedly applying the P-1 method
#
# Return the list of prime factors (in descending order)
func factorizeP1(n: BigInt, B1: int = 1000, B2: int = 100000): seq[BigInt] =
  var remaining = n
  while remaining > one:
    if isPrime(remaining):
      result.add(remaining)
      break

    # Remove small prime factors first
    var found = false
    for smallP in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
      let sp = smallP.initBigInt
      if remaining mod sp == zero:
        result.add(sp)
        remaining = remaining div sp
        found = true
        break

    if found:
      continue

    # Try the P-1 method
    let factor = pollardP1(remaining, B1 = B1, B2 = B2)

    if one < factor and factor < remaining:
      result.add(factor)
      remaining = remaining div factor
    else:
      result.add(remaining)  # left unfactored
      remaining = one
      break

  result.sort(order = Descending)

when isMainModule:
  let line = stdin.readLine
  let n = line.initBigInt
  echo "n = " & $n

  let factors = factorizeP1(n)
  let shown = "[" & join(factors.mapIt($it), ", ") & "]"
  echo "factors: " & shown
