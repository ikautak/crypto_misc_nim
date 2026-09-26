import strutils, os, math

type
  Vector = seq[float64]
  Matrix = seq[Vector]

proc dot(a, b: Vector): float64 =
  var sum = 0.0
  for i in 0 ..< a.len:
    sum += a[i] * b[i]
  return sum

proc sub(a, b: Vector): Vector =
  result = newSeq[float64](a.len)
  for i in 0 ..< a.len:
    result[i] = a[i] - b[i]
  return result

proc mul(a: Vector, scalar: float64): Vector =
  result = newSeq[float64](a.len)
  for i in 0 ..< a.len:
    result[i] = a[i] * scalar
  return result

proc gramSchmidt(b: Matrix): (Matrix, Matrix) =
  ## Returns (orthogonal basis u, mu coefficients)
  let n = b.len
  var u = newSeq[Vector](n)
  var mu = newSeq[Vector](n)
  for i in 0 ..< n:
    mu[i] = newSeq[float64](n)

  for i in 0 ..< n:
    u[i] = @b[i]
    for j in 0 ..< i:
      let u_j_dot_u_j = dot(u[j], u[j])
      if u_j_dot_u_j > 1e-12:
        mu[i][j] = dot(b[i], u[j]) / u_j_dot_u_j
        u[i] = sub(u[i], mul(u[j], mu[i][j]))
      else:
        mu[i][j] = 0.0
  
  return (u, mu)

proc lll(b: Matrix, delta: float64 = 0.75): Matrix =
  var B = newSeq[Vector](b.len)
  for i in 0 ..< b.len:
    B[i] = @b[i]
    
  let n = B.len
  var k = 1
  
  while k < n:
    var (u, _) = gramSchmidt(B)
    
    for j in countdown(k - 1, 0):
      let u_j_dot_u_j = dot(u[j], u[j])
      let mu_kj = if u_j_dot_u_j > 1e-12: dot(B[k], u[j]) / u_j_dot_u_j else: 0.0
      if abs(mu_kj) > 0.5:
        B[k] = sub(B[k], mul(B[j], round(mu_kj)))
        # Update GS basis because B[k] changed
        let (u_updated, _) = gramSchmidt(B)
        u = u_updated
        
    let (u_final, _) = gramSchmidt(B)
    let u_k = u_final[k]
    let u_kminus1 = u_final[k-1]
    let u_km1_dot_u_km1 = dot(u_kminus1, u_kminus1)
    let mu_k_kminus1 = if u_km1_dot_u_km1 > 1e-12: dot(B[k], u_kminus1) / u_km1_dot_u_km1 else: 0.0
    
    let lhs = dot(u_k, u_k)
    let rhs = (delta - mu_k_kminus1 * mu_k_kminus1) * u_km1_dot_u_km1
    
    if lhs >= rhs:
      k += 1
    else:
      let tmp = B[k]
      B[k] = B[k-1]
      B[k-1] = tmp
      k = max(1, k - 1)
      
  return B

proc readVector(line: string): Vector =
  # Split by comma or space
  var parts = line.split(',')
  if parts.len == 1:
    parts = line.split(' ')
    
  result = newSeq[float64]()
  for part in parts:
    let p = part.strip()
    if p != "":
      result.add(p.parseFloat())
  return result

proc main(filename: string) =
  if not fileExists(filename):
    echo "File not found: ", filename
    return

  var basis = newSeq[Vector]()
  for line in lines(filename):
    let trimmed = line.strip()
    if trimmed == "": continue
    basis.add(readVector(trimmed))
  
  if basis.len == 0:
    return

  let reduced = lll(basis)
  
  for v in reduced:
    var parts = newSeq[string]()
    for val in v:
      parts.add($int(round(val)))
    echo parts.join(", ")

when isMainModule:
  if paramCount() < 1:
    echo "Usage: lll <input_file>"
  else:
    main(paramStr(1))
