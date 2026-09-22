# Deterministic Miller-Rabin primality test.
#
# The witness set
#   { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 }
# is proven sufficient for every n < 2^64 (see "Deterministic variants of
# the Miller-Rabin primality test", Jim Sinclair / verified by others),
# It is not deterministic for n >= 2^64.

import bigints

const
  zero = 0.initBigInt
  one = 1.initBigInt
  two = 2.initBigInt

# One Miller-Rabin trial with witness a, where n - 1 = d * 2^r and d is odd.
# Returns true if n passes the trial (possibly prime).
func millerRabinTrial(a, d: BigInt, r: int, n: BigInt): bool =
  var x = powmod(a, d, n)
  if x == one or x == n - one:
    return true

  var i = 1
  while i < r:
    x = (x * x) mod n
    if x == n - one:
      return true
    if x == one: # non-trivial square root of 1 => n is composite
      return false
    inc i
  return false

func isPrime(n: BigInt): bool =
  if n < two:
    return false

  # Quick rejection / acceptance for small primes.
  let smallPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
  for p in smallPrimes:
    let pb = p.initBigInt
    if n == pb:
      return true
    if n mod pb == zero:
      return false

  # Write n - 1 = d * 2^r with d odd.
  var d = n - one
  var r = 0
  while d mod two == zero:
    d = d div two
    inc r

  # Witnesses sufficient for all n < 2^64 (deterministic).
  let witnesses = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
  for w in witnesses:
    let a = w.initBigInt mod n
    if a == zero:
      continue
    if not millerRabinTrial(a, d, r, n):
      return false
  return true

when isMainModule:
  let line = stdin.readLine()
  var n: BigInt
  try:
    n = line.initBigInt
  except ValueError:
    writeLine(stderr, "error: invalid input: " & line)
    quit(1)

  if n < zero:
    writeLine(stderr, "error: negative number: " & $n)
    quit(1)

  echo $n & ": " & (if isPrime(n): "prime" else: "composite")
