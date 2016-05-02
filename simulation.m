%This code simulate the phase reduction phenomenon in 3D settings assuming
%that the plane of focus is the center of the image
%function demo
    %object is placed at the central location
    clc;
    clear all;
    close all;
    %First, create a 3D block for simulation
    nx = 256;
    ny = 256;
    nz = 2000; %Make sure this number is large enough to obtain good sampling in z-domain
    uc = zeros(nx,ny,nz);%Total complex field
    %Object's dimension in microns
    ow = 30;%x
    oh = 30;%y
    ol = 25; %The thickness is 80 nm - z
    %Simulation area [um]
    sw = 80; %Dimension in the x
    sh = 80; %Dimension in the y
    sl = 80;
    %Simulation steps
    dx = sw/nx;
    dy = sh/ny;
    dz = sl/nz;
    %Wavelength in microns
    lambda = 0.574;beta0 = 2*pi/lambda;
    %Numberical aperture of the objectives
    NAc_arr = [0.0036,0.0072,0.014];
    m_phi_arr = zeros(nz,nx,length(NAc_arr));
    NAo = 0.35;%Estimated from Chris's data
    nbar = 1; %Free space refractive index
    n = 1.01; %Refractive index of quartz
    x_arr = linspace(-sw/2,sw/2,nx);
    y_arr = linspace(-sh/2,sh/2,ny);
    z_arr = linspace(-sl/2,sl/2,nz);
    u_arr = -nx/2:(nx/2-1);%Discrete frequency
    v_arr = -ny/2:(ny/2-1);
    w_arr = -nz/2:(nz/2-1);
    [u,v]=meshgrid(u_arr,v_arr);
    [x,y,z]=meshgrid(x_arr,y_arr,z_arr); %Generate the 3D coordinate
    xidx = intersect(find(x<ow/2),find(x>(-ow/2)));
    yidx = intersect(find(y<oh/2),find(y>(-oh/2)));
    zidx = intersect(find(z<ol/2),find(z>(-ol/2)));
    objidx = intersect(intersect(xidx,yidx),zidx);
    r = sqrt(x.^2+y.^2+z.^2);%Spherical object
    %objidx=find(r<10);
    kai = zeros(ny,nx,nz); %Note that the order is row - col and height or y, x, z 
    kai(objidx)=1;
    %Normalize to make sure that we have integration along z dimension will
    %gives us the object height
    clear xidx;clear yidx;clear zidx;
    kai = (n^2-nbar^2)*kai;
    figure(1);
    subplot(121);imagesc(x_arr,y_arr,kai(:,:,round(nz/2)),[0 max(kai(:))]);title('Kai - xy');xlabel('x [um]');ylabel('y [um]');colormap gray;colorbar
    subplot(122);xzview = squeeze(kai(round(ny/2),:,:));%Rotate the iamge so that it has
    xzview = xzview';
    imagesc(x_arr,z_arr,xzview,[0 max(kai(:))]);
    title('Kai - xz');xlabel('x [um]');ylabel('z [um]');colormap gray;colorbar
    
    kai = kai*dz;%if we have infinite number of samples kai*dz should be good here
    
    
    %Ukpz = zeros(nx,ny,nz); %Object field U(K_perp,z) at different planes
    Hkpz = zeros(nx,ny,nz); %Diffraction kernel at differnt z.
    %Next, compute the angular spectrum representation for each z values.
    %Note that the continuous wavevector k is related to the discrete
    %wavevector u by k = u/(M*deltax)
    %Compute the q map for the 2d map from the discrete map
    kx = 2*pi*u/nx/dx; %Continuous wavevector the unit is 1/micron
    ky = 2*pi*v/ny/dy;
    kp2 =(kx.^2+ky.^2);
    kx_arr = u_arr*2*pi/nx/dx;
    ky_arr = v_arr*2*pi/ny/dy;
    kz_arr = w_arr*2*pi/nz/dz;
    q2 =(nbar*beta0)^2-kp2; 
    mask = q2>=0;
    q=sqrt(q2).*mask;
    dq = (nbar*beta0-q).*mask;
  
    for NAcidx=1:length(NAc_arr)
        NAc = NAc_arr(NAcidx);
        %First, generate the convolution filter hi for the illumination
        kpthresh = NAc/nbar*beta0;
        kothresh = NAo/nbar*beta0; %Cutting frequency for the objective
        hof = kp2<kothresh^2;
        hif = exp(-kp2/2/kpthresh^2);
        hi = ifft2(ifftshift(hif));
        weighingcoef = sum(sum(hi)); %Normalize so that the filter has unit norm.
        hif = hif/weighingcoef;
        hi = ifft2(ifftshift(hif));
        figure(2);
        imagesc(kx_arr,ky_arr,(hif));title('hif');xlabel('kx (rad/um)');ylabel('ky (rad/um)');
        title('Condenser aperture profile k-domain');
        disp('Computing Kai(kpz)')
        tic;
        Kaikpz = fftshift(fftshift(fft(fft(kai,[],2),[],1),2),1);
        toc;
        %Next, compute the 3D Fourier transform 
        Kaikpkz =  fftshift(fft(Kaikpz,[],3),3);
        clear Ukpz;
        figure(3);
        subplot(121);imagesc(kx_arr,ky_arr,log10(abs(Kaikpkz(:,:,round(nz/2)))));title('Kai[kpkz] - xy view');xlabel('kx (rad/um)');ylabel('ky (rad/um)');
        subplot(122);xzview = squeeze(log10(abs(Kaikpkz(end/2,:,:))));xzview=xzview';%Rotate the iamge so that it has
        imagesc(kx_arr,kz_arr,xzview);title('Kai[kpkz] - xz view');xlabel('kx (rad/mu)');ylabel('kz (rad/um)');

        disp('Generating H[kpz]')
        %Compute the diffraction kernel for Hkpkz
        for zidx=1:nz
            Hkpz(:,:,zidx) = exp(1j*(dq)*z_arr(zidx)).*hof;     %Band limited due to the objective. The origin of z is now in the middle

        end
        %Normalize w.r.t to z. Make sure filtering in z doesn't get any
        %information on the depth
        Hkpkz = fftshift(fft(Hkpz,[],3),3); %
        phirpz = fftshift(real(beta0/2/nbar*ifftn(fftshift(Kaikpkz.*Hkpkz))),3); %Defocused phase, with no halo artifacts

         %Display the diffraction kernel
        figure(4);
        subplot(121);imagesc(kx_arr,ky_arr,log10(abs(Hkpkz(:,:,round(nz/2)))));title('Hkpkz - xy at view (at kz=0)');xlabel('kx (rad/um)');ylabel('ky (rad/um)');colorbar;
        subplot(122);xzview = squeeze(log10(abs(Hkpkz(end/2,:,:))));xzview=xzview';%Rotate the iamge so that it has
        imagesc(kx_arr,kz_arr,xzview);title('Hkpkz - xz view (at ky=0)');xlabel('kx (rad/um)');ylabel('kz (rad/um)');colorbar;

        z1 =10;%This is the plane of focus
        n1 = round(z1/dz);%Index difference from the pof to plane of interest
        z2 = -20;
        n2 = round(z2/dz);
        
        %Display the defocused phase
        figure(5);
        subplot(221);imagesc(x_arr,y_arr,phirpz(:,:,round(nz/2)));title('Phi (x,y,z=0)');xlabel('x [um]');ylabel('y [um]');colormap jet;colorbar;
        subplot(222);imagesc(x_arr,y_arr,phirpz(:,:,round(nz/2)-n1));title('Phi (x,y,z=5)');xlabel('x [um]');ylabel('y [um]');colormap jet;colorbar;
        subplot(223);imagesc(x_arr,y_arr,phirpz(:,:,round(nz/2)-n2));title('Phi (x,y,z=-20)');xlabel('x [um]');ylabel('y [um]');colormap jet;colorbar;
        subplot(224);imagesc(x_arr,z_arr,squeeze(phirpz(round(ny/2),:,:))');title('Phi (x,y=0,z)');xlabel('x [um]');ylabel('z [um]');colormap jet; colorbar;

         %Display the halo affected phase
        halo_phirpz = fftshift(real(beta0/2/nbar*ifftn(fftshift(Kaikpkz.*Hkpkz.*repmat((1-hif),[1 1 nz])))),3); %Phase under halo

        m_phi_arr(:,:,NAcidx) = squeeze(halo_phirpz(round(ny/2),:,:))'
        %Convert the phase into thickness in um
        Tdrpz = halo_phirpz/beta0/(n-nbar);
        for zidx=1:nz
            Tdrpz(:,:,zidx)=Tdrpz(:,:,zidx)-mean(mean(Tdrpz(1:20,round(nx/2),zidx)));
        end
        %Display the defocused phase
        figure(6);
        subplot(121);imagesc(x_arr,y_arr,halo_phirpz(:,:,round(nz/2)));title('Halo phi (x,y,z=0)');xlabel('x [um]');ylabel('y [um]');colorbar;
        subplot(122);imagesc(x_arr,z_arr,squeeze(halo_phirpz(round(ny/2),:,:))');title('Halo phi (x,y=0,z)');xlabel('x [um]');ylabel('z [um]');colorbar;



        figure(7);
        subplot(221);xzview = squeeze(Tdrpz(end/2,:,:));xzview = xzview';imagesc(x_arr,z_arr,xzview);title('Thickness - xz (nm)');xlabel('x (um)');ylabel('y (um)');colorbar;
        subplot(222);imagesc(x_arr,y_arr,Tdrpz(:,:,round(nz/2)-n2));title(strcat('Thickness - xy(nm), z = ',num2str(0),' um'));xlabel('x(um)');ylabel('y (um)');colorbar;
        subplot(223);imagesc(x_arr,y_arr,Tdrpz(:,:,round(nz/2)));title('Thickness - xy(nm), z = 0 um');xlabel('x(um)');ylabel('y (um)');colorbar;
        subplot(224);imagesc(x_arr,y_arr,Tdrpz(:,:,round(nz/2)+n2));title(strcat('Thickness - xy(nm), z = ',num2str(0),' um'));xlabel('x(um)');ylabel('y (um)');colorbar;
        figure(8);
        subplot(1,length(NAc_arr),NAcidx);
        imagesc(x_arr,z_arr,squeeze(m_phi_arr(:,:,NAcidx)));colormap jet;colorbar;xlabel('x [um]');ylabel('z [um]');
    end
  
     
    
    

