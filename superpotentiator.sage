# cyclic partial derivative of f with respect to the arrow a
def CPD(a, f):
  if len(f.terms()) > 1:
    return sum([CPD(a, term) for term in f.terms()])

  coefficient = f.leading_coefficient()
  term = f.support_of_term()

  result = 0
  for i in range(term.length()):
    cyclic = term[i:] * term[:i]

    if cyclic[0] == a.support_of_term():
      result = result + f.parent()(cyclic[1:])

  return coefficient * result


# compute the relations in the Jacobi algebra
def relations(w):
  return [CPD(w.parent(a), w) for a in w.parent().arrows()]


# algebra
Q = DiGraph({1: {1: ["x", "y", "z"]}})
kQ = Q.path_semigroup().algebra(QQ)
kQ.inject_variables()

# superpotential
w = x*y*z - x*z*y

print relations(w)
