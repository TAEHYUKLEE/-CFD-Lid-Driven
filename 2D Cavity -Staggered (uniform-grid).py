#!/usr/bin/env python
# coding: utf-8

# In[4]:


#2D Cavitry Staggered grids (uniform grids)
import numpy as np
from timeit import default_timer as timer
#import cupy as cp

#Global method
def TDMA(a,b,c,d,g,factor,nx,ny,nz,x,y,z):
    if x == 1: #x-direction TDMA
        for i in range(1,nx): #remove of a coefficient
            
            factor[i,:,:,:] = a[i,:,:,:]/b[i-1,:,:,:] #factor,a는 0 index는 버리는 것
            b[i,:,:,:] = b[i,:,:,:] - factor[i,:,:,:]*c[i-1,:,:,:]
            d[i,:,:,:] = d[i,:,:,:] - factor[i,:,:,:]*d[i-1,:,:,:]
            
        g[nx-1,:,:,:] = d[nx-1,:,:,:]/b[nx-1,:,:,:]
        
        for i in range(nx-2,(-1),-1):
            g[i,:,:,:] = (d[i,:,:,:]-c[i,:,:,:]*g[i+1,:,:,:])/b[i,:,:,:]
            
    elif y == 1: #y-direction TDMA
        for j in range(1,ny): #remove of a coefficient
            
            factor[:,j,:,:] = a[:,j,:,:]/b[:,j-1,:,:] 
            b[:,j,:,:] = b[:,j,:,:] - factor[:,j,:,:]*c[:,j-1,:,:]
            d[:,j,:,:] = d[:,j,:,:] - factor[:,j,:,:]*d[:,j-1,:,:]
            
        g[:,ny-1,:,:] = d[:,ny-1,:,:]/b[:,ny-1,:,:]
        
        for j in range(ny-2,(-1),-1):
            g[:,j,:,:] = (d[:,j,:,:]-c[:,j,:,:]*g[:,j+1,:,:])/b[:,j,:,:]
    
    elif z == 1: #y-direction TDMA
        for k in range(1,nz): #remove of a coefficient
            
            factor[:,:,k,:] = a[:,:,k,:]/b[:,:,k-1,:] 
            b[:,:,k,:] = b[:,:,k,:] - factor[:,:,k,:]*c[:,:,k-1,:]
            d[:,:,k,:] = d[:,:,k,:] - factor[:,:,k,:]*d[:,:,k-1,:]
            
        g[:,:,nz-1,:] = d[:,:,nz-1,:]/b[:,:,nz-1,:]
        
        for k in range(nz-2,(-1),-1):
            g[:,:,k,:] = (d[:,:,k,:]-c[:,:,k,:]*g[:,:,k+1,:])/b[:,:,k,:]
    
    return g

def Poisson(p_n,p_o,u_mid,v_mid,dt,nx,ny,dx,dy,epsil_p):
    
    steps = int(0)
    
    while(1):
        
        sum_error = np.zeros([1],dtype = np.float64)
        
        if steps ==0:
            p_o[0,1:ny,:,:] = p_o[1,1:ny,:,:]; p_o[nx,1:ny,:,:] = p_o[nx-1,1:ny,:,:]
            p_o[1:nx,0,:,:] = p_o[1:nx,1,:,:]; p_o[1:nx,ny,:,:] = p_o[1:nx,ny-1,:,:]  
        
        p_n[1:nx,1:ny,:,:] = (p_o[2:nx+1,1:ny,:,:] + p_o[0:nx-1,1:ny,:,:]+p_o[1:nx,2:ny+1,:,:] + p_o[1:nx,0:ny-1,:,:] - dx**2/dt*                             ((u_mid[2:nx+1,1:ny,:]-u_mid[1:nx,1:ny,:])/dx +(v_mid[1:nx,2:ny+1,:]-v_mid[1:nx,1:ny,:])/dy))/4.0
        
        #Neumman B.C
        p_n[0,1:ny,:,:] = p_n[1,1:ny,:,:]; p_n[nx,1:ny,:,:] = p_n[nx-1,1:ny,:,:]
        p_n[1:nx,0,:,:] = p_n[1:nx,1,:,:]; p_n[1:nx,ny,:,:] = p_n[1:nx,ny-1,:,:]
        
        sum_error[:] = np.sum(abs(p_n[1:nx,1:ny] - p_o[1:nx,1:ny]))/((nx-1)*(ny-1))
        
        p_o[:]= p_n[:]
        
        if (sum_error<epsil_p):
            break
            
        steps = steps+1
        
    return p_o, steps

def main():
    
    #grids shape
    nx = 129; ny = 129
    
    #---------------------- Define variables ---------------------#
    #Dfine 2-D velocity
    u_o = np.zeros([nx+1,ny+1,1,1], dtype = np.float64);v_o = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    u_n = np.zeros([nx+1,ny+1,1,1], dtype = np.float64);v_n = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    u_c = np.zeros([nx,ny,1,2], dtype = np.float64)
    stream = np.zeros([nx,ny,1,1], dtype = np.float64)
    vorticity = np.zeros([nx-1,ny-1,1,1],dtype = np.float64)
    p_o = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    p_n = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    
    #Coefficient
    #x-direction coefficient
    a_x = np.zeros([nx+1,1,1,1],dtype = np.float64); b_x = np.zeros([nx+1,1,1,1],dtype = np.float64)
    c_x = np.zeros([nx+1,1,1,1],dtype = np.float64); d_u = np.zeros([nx-2,ny-1,1,1],dtype = np.float64)
    factor_x = np.zeros([nx+1,1,1,1],dtype = np.float64)

    #y-direction coefficient
    a_y = np.zeros([1,ny+1,1,1],dtype = np.float64); b_y = np.zeros([1,ny+1,1,1],dtype = np.float64)
    c_y = np.zeros([1,ny+1,1,1],dtype = np.float64); d_v = np.zeros([nx-1,ny-2,1,1],dtype = np.float64)
    factor_y = np.zeros([1,ny+1,1,1],dtype = np.float64)
    
    #intermediate factor for TDMA
    g = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    h = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    u_mid = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    v_mid = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    
    #gradient of velocity
    dudx = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); dvdx = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    dudy = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); dvdy = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    dudx_2 = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); dvdx_2 = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    dudy_2 = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); dvdy_2 = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    
    #Non-linear
    H_uo = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); H_vo = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    H_up = np.zeros([nx+1,ny+1,1,1], dtype = np.float64); H_vp = np.zeros([nx+1,ny+1,1,1], dtype = np.float64)
    #-------------------------------------------------------------#
    
    # Simulation parameter
    Re = 50000
    dt = 0.0025
    t_real =  int(0)
    epsil_p = 10**(-8)
    epsil_u = 10**(-8)
    
    #simulation factor
    x_len = 1.0; y_len = 1.0
    
    x = np.zeros([nx,1,1,1], dtype = np.float64)
    y = np.zeros([1,ny,1,1], dtype = np.float64)
    
    dx = x_len/(nx-1); dy = y_len/(ny-1) #delta panel 칸의 갯수로 나눠줘야함
    x0 = 0.0; y0 = 0.0
    xf = 1.0; yf = 1.0
    
    alpha = 1.0
    video_scene=0
    
    #uniform grids
    x[:,0,0,0] = np.linspace(x0,xf,nx)
    y[0,:,0,0] = np.linspace(y0,yf,ny)
    
    mu = dt/(2*Re*dx**2)
    
    start = timer()
    
    while(1):
        
        #update non-linear term
        H_up[2:nx,1:ny] = H_uo[2:nx,1:ny]
        H_vp[1:nx,2:ny] = H_vo[1:nx,2:ny]
        
        #Boundary condition
        u_o[1:nx,ny,:,:] = (2 - u_o[1:nx,ny-1,:,:]); u_o[1:nx,0,:,:] = (2-u_o[1:nx,1,:,:])
        v_o[nx,1:ny,:,:] = (2- v_o[nx-1,1:ny,:,:]); v_o[0,1:ny,:,:] = (2- v_o[1,1:ny,:,:])
    
        dudx[2:nx,1:ny,:,:] = 0.5*((u_o[3:nx+1,1:ny,:,:] - u_o[2:nx,1:ny,:,:])/dx + (u_o[2:nx,1:ny,:,:] - u_o[1:nx-1,1:ny,:,:])/dx)
        dudy[1:nx,1:ny,:,:] = 0.5*((u_o[1:nx,2:ny+1,:,:] - u_o[1:nx,1:ny,:,:])/dy + (u_o[1:nx,1:ny,:,:] - u_o[1:nx,0:ny-1,:,:])/dy)
        #dudx 범위는 2:nx, 1:ny로 정해진다
        dvdx[1:nx,1:ny,:,:] = 0.5*((v_o[2:nx+1,1:ny,:,:] - v_o[1:nx,1:ny,:,:])/dx + (v_o[1:nx,1:ny,:,:] - v_o[0:nx-1,1:ny,:,:])/dx)
        dvdy[1:nx,2:ny,:,:] = 0.5*((v_o[1:nx,3:ny+1,:,:] - v_o[1:nx,2:ny,:,:])/dy + (v_o[1:nx,2:ny,:,:] - v_o[1:nx,1:ny-1,:,:])/dy)
        
        H_uo[2:nx,1:ny] = u_o[2:nx,1:ny]*dudx[2:nx,1:ny]+0.5*(0.5*(v_o[1:nx-1,1:ny] + v_o[1:nx-1,2:ny+1])+0.5*(v_o[2:nx,1:ny] + v_o[2:nx,2:ny+1]))*dudy[2:nx,1:ny]
        H_vo[1:nx,2:ny] = 0.5*(0.5*(u_o[1:nx,1:ny-1]+u_o[1:nx,2:ny]) + 0.5*(u_o[2:nx+1,1:ny-1]+u_o[2:nx+1,2:ny]))*dvdx[1:nx,2:ny] + v_o[1:nx,2:ny]*dvdy[1:nx,2:ny]
        
        dudx_2[2:nx,1:ny] = ((u_o[3:nx+1,1:ny] - u_o[2:nx,1:ny])/dx - (u_o[2:nx,1:ny] - u_o[1:nx-1,1:ny])/dx)/dx
        dudy_2[1:nx,1:ny] = ((u_o[1:nx,2:ny+1] - u_o[1:nx,1:ny])/dy - (u_o[1:nx,1:ny] - u_o[1:nx,0:ny-1])/dy)/dy
        dvdx_2[1:nx,1:ny] = ((v_o[2:nx+1,1:ny] - v_o[1:nx,1:ny])/dx - (v_o[1:nx,1:ny] - v_o[0:nx-1,1:ny])/dx)/dx
        dvdy_2[1:nx,2:ny] = ((v_o[1:nx,3:ny+1] - v_o[1:nx,2:ny])/dy - (v_o[1:nx,2:ny] - v_o[1:nx,1:ny-1])/dy)/dy
        #Debugging 여기에 0.5가 곱해져 있었음
        
        if t_real ==0: H_up[2:nx,1:ny] = H_uo[2:nx,1:ny]; H_vp[1:nx,2:ny] = H_vo[1:nx,2:ny] #only first step
            
        #d coefficient (Boundary condition은 따로 처리해줘야 한다)
        d_u[:]= dt/(Re)*(dudx_2[2:nx,1:ny] + dudy_2[2:nx,1:ny]) - dt/2*(3*H_uo[2:nx,1:ny] - H_up[2:nx,1:ny]) 
        d_v[:]= dt/(Re)*(dvdx_2[1:nx,2:ny] + dvdy_2[1:nx,2:ny]) - dt/2*(3*H_vo[1:nx,2:ny] - H_vp[1:nx,2:ny])
        
        #------------- u TDMA  ---------------#
        #x-direction TDMA for u velocity
        a_x[2:nx] = -mu; b_x[2:nx] = 1+2*mu; c_x[2:nx] = -mu
        
        x_d=1; y_d=0; z_d=0
        g[2:nx,1:ny,:,:] = TDMA(a_x[2:nx],b_x[2:nx],c_x[2:nx],d_u,g[2:nx,1:ny],factor_x[2:nx],nx-2,ny-1,1,x_d,y_d,z_d)

        #y-direction TDMA for u velocity
        a_y[:,1:ny] = -mu; b_y[:,1:ny] = 1+2*mu; c_y[:,1:ny] = -mu
        b_y[:,1] = (1+3*mu); b_y[:,ny-1] = (1+3*mu)
        
        x_d=0; y_d=1; z_d=0
        u_mid[2:nx,1:ny,:,:] = TDMA(a_y[:,1:ny],b_y[:,1:ny],c_y[:,1:ny],g[2:nx,1:ny],h[2:nx,1:ny],factor_y[:,1:ny],nx-2,ny-1,1,x_d,y_d,z_d)
        
        #------------ v TDMA -----------------#
        #x-direction TDMA for v velocity
        a_x[1:nx] = -mu; b_x[1:nx] = 1+2*mu; c_x[1:nx] = -mu
        b_x[1,:] = (1+3*mu); b_x[nx-1,:] = (1+3*mu)
        
        x_d=1; y_d=0; z_d=0
        g[1:nx,2:ny,:,:] = TDMA(a_x[1:nx],b_x[1:nx],c_x[1:nx],d_v,g[1:nx,2:ny],factor_x[1:nx],nx-1,ny-2,1,x_d,y_d,z_d)

        #y-direction TDMA for v velocity
        a_y[:,2:ny] = -mu; b_y[:,2:ny] = 1+2*mu; c_y[:,2:ny] = -mu
        
        x_d=0; y_d=1; z_d=0
        v_mid[1:nx,2:ny,:,:] = TDMA(a_y[:,2:ny],b_y[:,2:ny],c_y[:,2:ny],g[1:nx,2:ny],h[1:nx,2:ny],factor_y[:,2:ny],nx-1,ny-2,1,x_d,y_d,z_d)
        
        
        #update intermeidate velocity
        u_mid[2:nx,1:ny] = u_o[2:nx,1:ny] + u_mid[2:nx,1:ny]
        v_mid[1:nx,2:ny] = v_o[1:nx,2:ny] + v_mid[1:nx,2:ny]
        
        
        #-------------- Poisson equation ----------------#
        p_o, steps = Poisson(p_n,p_o,u_mid,v_mid,dt,nx,ny,dx,dy,epsil_p)
        #Poisson equation에서 update된 p_o를 받아와야 하는데 p_n을 받아와서 터진적이 있음
        
        #Correction pressure term
        u_n[2:nx,1:ny,:,:] = u_mid[2:nx,1:ny,:,:] - dt*(p_o[2:nx,1:ny] - p_o[1:nx-1,1:ny])/dx #u_update
        v_n[1:nx,2:ny,:,:] = v_mid[1:nx,2:ny,:,:] - dt*(p_o[1:nx,2:ny] - p_o[1:nx,1:ny-1])/dy #v_update
        
        print("Residual")
        print(np.sum(abs(u_n[1:nx,1:ny] - u_o[1:nx,1:ny]) + abs(v_n[1:nx,1:ny] - v_o[1:nx,1:ny])/((nx-1)*(ny-1))))
        
        Residual = (np.sum(abs(u_n[1:nx,1:ny] - u_o[1:nx,1:ny]) + abs(v_n[1:nx,1:ny] - v_o[1:nx,1:ny])/((nx-1)*(ny-1))))
        
        if (Residual < epsil_u):
            break
        

        u_o[:] = u_n[:]; v_o[:] = v_n[:]
        
        
        #print contour stream line & velocity field
        '''
        if t_real% 20==0:
            video_plt = open('2D u_final_video.plt', 'a')

            if t_real==0: video_plt.write('VARIABLES = "x","y","u","v"\n')
            video_plt.write('Zone T="HIT%d"\n'%0)
            video_plt.write('I=%d J=%d\n' %(nx-1, ny-1))
    
            for j in range(1, ny):
                for i in range(1, nx):
                    video_plt.write('%.5f %.5f %.5f %.5f\n' %(i,j, u_o[i,j,0,0], v_o[i,j,0,0]))
    
            video_plt.close()
        '''

            #Boundary condition
            #u_n[1:nx,ny,:,:] = (2 - u_n[1:nx,ny-1,:,:]); u_n[1:nx,0,:,:] = -u_n[1:nx,1,:,:]
            #v_n[nx,1:ny,:,:] = - v_n[nx-1,1:ny,:,:]; v_n[0,1:ny,:,:] = - v_n[1,1:ny,:,:]
            
        if t_real%50==0:
            
            #Staggered grids to Collocated grids
            u_c[0:nx,0:ny,:,0:1] = (u_n[1:nx+1,0:ny] + u_n[1:nx+1,1:ny+1])/2
            u_c[0:nx,0:ny,:,1:2] = (v_n[0:nx,1:ny+1] + v_n[1:nx+1,1:ny+1])/2   
            
            #Get vorticity [d(v_y)/dx - d(v_x)/dy]
        
            vorticity[:] = ((u_c[1:nx,1:ny,:,1:2]-u_c[0:nx-1,1:ny,:,1:2])/dx + (u_c[1:nx,0:ny-1,:,1:2]-u_c[0:nx-1,0:ny-1,:,1:2])/dx)/2-            ((u_c[0:nx-1,1:ny,:,0:1] - u_c[0:nx-1,0:ny-1,:,0:1])/dy + (u_c[1:nx,1:ny,:,0:1] - u_c[1:nx,0:ny-1,:,0:1])/dy)/2
            
            '''
            for i in range(nx):
        
                sum_stream = np.zeros([1], dtype = np.float64)
        
                for j in range(ny):
                    sum_stream[:] = sum_stream[:] + u_c[i,j,:,0]*dy
                    stream[i,j] = sum_stream[:]
                    
            stream_video = open('2D stream-line.plt', 'a')

            if t_real==0: stream_video.write('VARIABLES = "x","y","stream"\n')
            stream_video.write('Zone T="HIT%d"\n'%0)
            stream_video.write('I=%d J=%d\n' %(nx, ny))
    
            for j in range(ny):
                for i in range(nx):
                    stream_video.write('%.5f %.5f %.5f\n' %(i,j, stream[i,j,0,0]))
                    
            stream_video.close()     
            '''
            
            vorticity_video = open('2D vorticity.plt', 'a')
            
            if t_real==0: vorticity_video.write('VARIABLES = "x","y","vorticity"\n')
            vorticity_video.write('Zone T="HIT%d"\n'%0)
            vorticity_video.write('I=%d J=%d\n' %(nx-1, ny-1))
    
            for j in range(ny-1):
                for i in range(nx-1):
                    vorticity_video.write('%.5f %.5f %.5f\n' %(i,j, vorticity[i,j,0,0]))
                    
            vorticity_video.close()  
            
            video_scene += 1 
            
            if video_scene == 410: break
            print(video_scene)
            
        t_real += 1               
            
'''        
        
     
    
    duration = timer() - start 
        
    #Boundary condition
    u_n[1:nx,ny,:,:] = (2 - u_n[1:nx,ny-1,:,:]); u_n[1:nx,0,:,:] = -u_n[1:nx,1,:,:]
    v_n[nx,1:ny,:,:] = - v_n[nx-1,1:ny,:,:]; v_n[0,1:ny,:,:] = - v_n[1,1:ny,:,:]
    
    #Staggered grids to Collocated grids
    u_c[0:nx,0:ny,:,0:1] = (u_n[1:nx+1,0:ny] + u_n[1:nx+1,1:ny+1])/2
    u_c[0:nx,0:ny,:,1:2] = (v_n[0:nx,1:ny+1] + v_n[1:nx+1,1:ny+1])/2    
    
    #stream line computation
    for i in range(nx):
        
        sum_stream = np.zeros([1], dtype = np.float64)
        
        for j in range(ny):
            sum_stream[:] = sum_stream[:] + u_c[i,j,:,0]*dy
            stream[i,j] = sum_stream[:]
            
    
    #print contour stream line & velocity field
    contour = open('2D u_final.plt', 'w')

    contour.write('VARIABLES = "x","y","u","v"\n')
    contour.write('Zone T="HIT%d"\n'%0)
    contour.write('I=%d J=%d\n' %(nx-1, ny-1))
    
    for j in range(1, ny):
        for i in range(1, nx):
            contour.write('%.5f %.5f %.5f %.5f\n' %(i,j, u_o[i,j,0,0], v_o[i,j,0,0]))
    
    contour.close()
    
    print("end")
    print(duration)
    
    #print stream line & Collocated grids
    
    contour = open('2D stream-line.plt', 'w')

    contour.write('VARIABLES = "x","y","stream"\n')
    contour.write('Zone T="HIT%d"\n'%0)
    contour.write('I=%d J=%d\n' %(nx, ny))
    
    for j in range(ny):
        for i in range(nx):
            contour.write('%.5f %.5f %.5f\n' %(i,j, stream[i,j,0,0]))
    
    contour.close()    
                
#vorticity (x-direction) d(v_z)/dy - d(v_y)/dz
vorticity_gpu[:,:,:,0] = (up_storage_gpu[1:nyp,1:nzp,1:nxp,2] - up_storage_gpu[0:nyp-1,1:nzp,1:nxp,2])/dy_gpu[0:nyp-1, :, :,0] \
                                        - (up_storage_gpu[1:nyp,1:nzp,1:nxp,1] - up_storage_gpu[1:nyp,0:nzp-1,1:nxp,1])/dz

#vorticity (y-direction) d(v_x)/dz - d(v_z)/dx
vorticity_gpu[:,:,:,1] = (up_storage_gpu[1:nyp,1:nzp,1:nxp,0] - up_storage_gpu[1:nyp, 0:nzp-1,1:nxp,0])/dz \
                                - (up_storage_gpu[1:nyp,1:nzp,1:nxp,2] - up_storage_gpu[1:nyp,1:nzp,0:nxp-1,2])/dx

#vorticity (z-direction) d(v_y)/dx - d(v_x)/dy
vorticity_gpu[:,:,:,2] = (up_storage_gpu[1:nyp,1:nzp,1:nxp,1] - up_storage_gpu[1:nyp,1:nzp,0:nxp-1,1])/dx \
                    - (up_storage_gpu[1:nyp,1:nzp,1:nxp,0] - up_storage_gpu[0:nyp-1,1:nzp,1:nxp,0])/dy_gpu[0:nyp-1, :, :,0]        
        
'''    

            

        
if __name__ == "__main__":
        
    main()
        
    print("All process is finished!")


# In[ ]:




