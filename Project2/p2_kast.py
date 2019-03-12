for j in range(N):
    dW_dx = lambda W_val,x: 1/x**2*(np.exp(2*W_eq[j] - W_val) - np.exp(W_val))  #overflow
    k1 = dW_dx(W[j,ind], x_arr[j])
    k2 = dW_dx(W[j,ind] + k1*h/2., x_arr[j] + h/2.)
    k3 = dW_dx(W[j,ind] + k2*h/2., x_arr[j] + h/2.)
    k4 = dW_dx(W[j,ind] + k3*h, x_arr[j] + h)
    W[j+1,ind] = W[j,ind] + h*(k1 + 2*k2 + 2*k3 + k4)/6. #invalid value encountered in double_scalars
    #W[j+1,ind] = W[j,ind] + h*dW_dx(W[j,ind] + h*dW_dx(W[j,ind],x_arr[j])/2., x_arr[j] + h/2.)
    print W[j,ind], yeq[j], W_eq[j]
#print W[:,ind]
W_eq_arr[:,ind] = W_eq
ind += 1
