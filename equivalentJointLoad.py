def equivalentJointLoad(w):
    choice=int(input("Enter 1: UDL\t2: Point Load at Mid Length\t3: Point load applied at arbitory point on the beam\t0: If no such load applied: "))
    
    #For uniformly distributed load.
    if choice==1:
        #w=float(input("Enter the value of load being applied: "))
        l=float(input("Enter the length of beam applied on: "))
        moment1= ((w*l**2)/12)
        moment2= ((w*l**2)/12)
        force1 = ((w*l)/2)
        force2 = ((w*l)/2)
        
        eq=moment1,moment2,force1,force2
        
        
    elif choice==2:
        #w=float(input("Enter the value of load being applied: "))
        l=float(input("Enter the length of beam applied on: "))
        moment1= ((w*l)/8)
        moment2= ((w*l)/8)
        force1 = ((w)/2)
        force2 = ((w)/2)
        
        eq=moment1,moment2,force1,force2
    
    elif choice==0:
        moment1= 0
        moment2= 0
        force1 = 0
        force2 = 0
        eq=moment1,moment2,force1,force2
    
    elif choice==3:
        w=float(input("Enter the value of load being applied: "))
        a=float(input("Enter the distance of load from left support: "))
        b=float(input("Enter the distance of load from right support: "))
        l=a+b
        moment1= ((w*a*b**2)/l**2)
        moment2= ((w*a**2*b)/l**2)
        force1 = ((w*l)/2)
        force2 = ((w*l)/2)
    else:
        print("Enter the equivalent joint load properly.")
        
    return eq