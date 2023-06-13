# X = [1000, 1500, 2000, 3000, 4000, 5000, 6000]
# Y = [ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500]
# YP =[ 30,  60,  90, 120, 150, 200, 250, 300, 400, 500]

def mass_points():

    mpoints=[]
    mYprime_l=30
    mYprime_r_1=250
    mYprime_r_2=500
    mYprime_step1=30
    mYprime_step2=50
    mYprime_step3=100
    mY_l=30
    mY_step1=30
    mY_step2=50
    mY_step3=100
    mX_l=1000
    mX_r=6000
    mX_step1=500
    mX_step2=1000

    mX_tmp=mX_l
    mX_step_tmp=mX_step1
    mYprime_r_tmp=mYprime_r_1
    while mX_tmp <= mX_r:
            mYprime_tmp=mYprime_l
            mYprime_step_tmp=mYprime_step1

            while mYprime_tmp <= mYprime_r_tmp:
                mY_tmp=mY_l
                mY_step_tmp=mY_step1
                
                while mY_tmp <= mYprime_tmp:

                    mpoints.append([mX_tmp,mY_tmp,mYprime_tmp])

                    if mY_tmp == 150:
                        mY_step_tmp=mY_step2
                    elif mY_tmp == 300:
                        mY_step_tmp=mY_step3

                    mY_tmp+=mY_step_tmp
            
                if mYprime_tmp == 150:
                    mYprime_step_tmp=mYprime_step2
                elif mYprime_tmp == 300:
                    mYprime_step_tmp=mYprime_step3

                mYprime_tmp+=mYprime_step_tmp

            if mX_tmp == 2000:
                mX_step_tmp=mX_step2
                mYprime_r_tmp=mYprime_r_2

            mX_tmp+=mX_step_tmp
    
    return map(tuple,mpoints)