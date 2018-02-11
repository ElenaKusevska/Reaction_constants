program partition
implicit none

!---------------------------------------
! calculates the rate constant
! for the reaction H + H2 -> H-H-H  
!---------------------------------------

integer, parameter :: dp = SELECTED_REAL_KIND(15)
integer :: i !counter

!------------
!constants:
!------------

real(dp), parameter :: boltzmann = 1.380648813D-23 ![m^(2)kg/(s^(2)K)]
real(dp), parameter :: h = 6.626069573D-34 ![m^(2)kg/s]
real(dp), parameter :: pi = 3.141592654
real(dp), parameter :: c = 299792458  !speed of light - [m/s]
real(dp), parameter :: Ea = 9.7 ![kcal/mol]
real(dp), parameter :: kcaltoJ = 4184 !conversion factor
real(dp), parameter :: angstrom = 1D-10 !convert angstrom to meters
real(dp), parameter :: mp = 1.672621777D-27 !mass of proton - [kg]
real(dp), parameter :: me = 9.10938291D-31 !mass of electron - [kg]
real(dp), parameter :: symm = 2 !symmetry number.
real(dp), parameter :: R = 8.314462175 !gas constant [J/(Kmol)]
real(dp), parameter :: N = 6.0221412927D+23!avogadro's number - [mole^(-1)]

!---------------
!variables:
!---------------

real, parameter :: V = 1 !volume -> [dm^3]
real(dp) :: E, vim, elec !Energy to cross the barrier, imaginary
                              !frequency and electronic p.f.
real(dp), dimension(701) :: k, Q, T, trans, vib, rot! rate constant, 
                                       !temperature and partition 
                                       !function for the reaction
real(dp), parameter :: eH = 2, eH2 = 1, eHHH = 2 !electronic partition 
                                                      !functions 
real(dp), dimension(701) :: dG, dH, dS, KC, wigner, kcorr 
                                 !Gibbs free energy, enthalpy and 
                                 !entropy of the transition state
                                 !reaction, the equilibrium constant and 
                                 ! the wigner correction.                  

!Hydrogen:
real(dp) :: mH !mass of hydrogen
real(dp), dimension(701) :: qtrans

!Hydrogen molecule:
real(dp), dimension(701) :: vibH2, rotH2, transH2
real(dp) :: f !vibrational frequency [s^(-1)].
real(dp) :: mH2, rH2 !mass of hydrogen atom, mass 
                                 !and radius of hydrogen molecule
real(dp) :: IH2 ! angular momentum
real(dp) :: ZPEH2 !ZP correction

!H-H-H:
real(dp), dimension(701) :: HHHrot, HHHvib, HHHtrans, QHHH
real(dp) :: IHHH, rHHH, mHHH !rotational momentum, distances and mass for 
                                             !the transition state
real(dp) :: vHHH1, vHHH2 ! vibrational frequencies, ali toa vtoroto 
                                       !moze da bide greska
real(dp) :: ZPEHHH !zero point correction


!---------------------------------
! Define some values first:
!---------------------------------

vim = 850*c *100
elec = eHHH/(eH2*eH)
mH = mp + me ![kg] - mass of hydrogen
f = 4401 * c * 100  
mH2 = (2.)*(mp + me) ![kg]
rH2 = 0.741*angstrom !distance between hydrogen atoms - [m]
IH2 = ((mH*mH)/(mH+mH))*(rH2**(2)) !moment of inertia of hydrogen molecule
ZPEH2 = (h*f*N)/2


!moment of inertia and mass of transition state:
rHHH = 0.929*angstrom
IHHH = (2.)*mH*(rHHH**(2.))
!IDHH = mD*(rDHH**(2)) + mH*(rDHH**(2)) - ((mD*rDHH-mH*rDHH)**(2))/(mD+mH+mH)
mHHH = mH + mH + mH
vHHH1 = (1780.)*c*(100.)
vHHH2 = (861.)*c*(100.)
ZPEHHH = Ea*kcaltoJ + N*(0.5)*(h*vHHH1 + (2.)*h*vHHH2) 

E = (ZPEHHH-ZPEH2)

open (unit=1, file='results_H2.txt', status='new', action='write')
write(1,*) 'T', '	', 'k', '	', 'Kc', '	', 'G', '	', 'H', '	', 'S', '	',&
'wigner', '	', 'kcorr'

!-----------------------------------------
! calculate the partition functions:
!-----------------------------------------

do i = 1, 701 !701  
   T(i) = i + 299

   !-------------------------------------------------------------Deuterium:
   !-----------------------------------------------------------------------

   qtrans(i) = ((((2.)*pi*mH*boltzmann*T(i))**(.3/.2))/(h**(3.)))*V

   !-----------------------------------------------------Hydrogen molecule:
   !-----------------------------------------------------------------------

   vibH2(i) = (1.) / ((1.) - exp(-h*f/(boltzmann*T(i))))
   rotH2(i) = ((8.)*(pi**(2.))*IH2*boltzmann*T(i))/h
   transH2(i) = (((2.)*pi*mH2*boltzmann*T(i))**(3./2.))/(h**(3.))*V

   !-----------------------------------------------------------------D-H-H:
   !-----------------------------------------------------------------------

   HHHrot(i) = ((8.)*(pi**(2.))*IHHH*boltzmann*T(i))/h
   HHHvib(i) = ((1.) / ((1.) - exp(-h*vHHH2/(boltzmann*T(i))))**(2.))*((1.) / ((1.) - exp(-h*vHHH1/(boltzmann*T(i)))))

   HHHtrans(i) = (((2.)*pi*mHHH*boltzmann*T(i))**(3./2.))/(h**(3.))*V

   !------------------------------------------------------------value of k:
   !-----------------------------------------------------------------------

   rot(i) = HHHrot(i)/rotH2(i)
   trans(i) = HHHtrans(i)/(qtrans(i)*transH2(i))
   vib(i) = HHHvib(i)/vibH2(i)   
   Q(i) = rot(i)*vib(i)*trans(i)*elec
   KC(i) = Q(i)*N*(1000.)*exp(-(E/(R*T(i)))) !L/mol
   k(i) = symm*((boltzmann*T(i))/h)*Q(i)*exp(-(E/(R*T(i))))*N*1000
   dG(i) = -R*T(i)*LOG(KC(i))
   dH(i) = E - (2.)*R*T(i)
   dS(i) = (dH(i)-dG(i))/T(i)
   wigner(i) = (1.) + ((1./24.))*((h*vim)/(boltzmann*T(i)))**(2.)
   kcorr(i) = k(i)*wigner(i)

   write(1,*) T(i), k(i), Kc(i), dG(i), dH(i), dS(i), wigner(i), kcorr(i)
end do

end program partition
