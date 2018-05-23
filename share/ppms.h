
#define   I_DOMAIN        I1-1:IN+1

#define range(a,b) I1-1+a:IN-1+a, J1-1+b:JN-1+b

#define   BOUNDS           I1:IN

#define  PADDED_BOUNDS     I1-5:IN+5

#define  PAD_1_BOUNDS     I1-1:IN+1

! ASSUMING CONVENTION
!    
!          ---- F ---- Q ---- F ---- Q ---
!               i      i     i+1    i+1


#define  F_1_RT_OF_Q  I1+1:IN+1
#define  F_1_LF_OF_Q  I1:IN

#define  F_2_RT_OF_Q  I1+2:IN+2
#define  F_2_LF_OF_Q  I1-1:IN-1

#define  F_3_RT_OF_Q  I1+3:IN+3
#define  F_3_LF_OF_Q  I1-2:IN-2



#define  Q_1_RT_OF_F  I1:IN
#define  Q_1_LF_OF_F  I1-1:IN-1

#define  Q_2_RT_OF_F  I1+1:IN+1
#define  Q_2_LF_OF_F  I1-2:IN-2

#define  Q_3_RT_OF_F  I1+2:IN+3
#define  Q_3_LF_OF_F  I1-3:IN-3




#define  X_1_RT_OF_X  I1+1:IN+1
#define  X_1_LF_OF_X  I1-1:IN-1

#define  X_2_RT_OF_X  I1+2:IN+2
#define  X_2_LF_OF_X  I1-2:IN-2






