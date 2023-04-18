KPL/MK
 
      Meta-kernel for the ``Mars MEX HRSC camera data pre handle '' task in the file
      Hands On Lesson.
 
      The names and contents of the kernels referenced by this
      meta-kernel are as follows:
 
      File name                    			Contents
      ---------------------------  			------------------------------------
      NAIF0012.TLS                 			Generic LSK
	  MEX_210907_STEP.TSC		   			MEX    SCLK
	  ATNM_MEASURED_040101_050101_V03.BC	MEX CK
	  PCK00010.TPC				   			planet PCK
      DE405.BSP          		   			Solar System Ephemeris
	  ORMM__031222180906_00052.BSP 			MEX Ephemeris
      MEX_V16.TF               	   			MEX instrument FK
 
      \begindata
 
         KERNELS_TO_LOAD = ( 'D:\\Code\\Mars\\data\\kernel\\NAIF0012.TLS'
                             'D:\\Code\\Mars\\data\\kernel\\MEX_210907_STEP.TSC'
                             'D:\\Code\\Mars\\data\\kernel\\MEX_V16.TF'
                             'D:\\Code\\Mars\\data\\kernel\\ATNM_MEASURED_040101_050101_V03.BC'
							 'D:\\Code\\Mars\\data\\kernel\\ATNM_MEASURED_2013_V04.BC'
							 'D:\\Code\\Mars\\data\\kernel\\ATNM_MEASURED_2017_V01.BC'
							 'D:\\Code\\Mars\\data\\kernel\\ATNM_MEASURED_180101_181231_V01.BC'
							 'D:\\Code\\Mars\\data\\kernel\\ATNM_MEASURED_2019_V02.BC'
                             'D:\\Code\\Mars\\data\\kernel\\DE405.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\PCK00010.TPC'
                             'D:\\Code\\Mars\\data\\kernel\\ORMM__031222180906_00052.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM__131201000000_01033.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM_T19_170101000000_01318.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM_T19_170201000000_01326.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM_T19_170301000000_01336.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM_T19_170401000000_01344.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM__181101000000_01484.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM__181201000000_01490.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\ORMM_T19_190101000000_01498.BSP'
							 'D:\\Code\\Mars\\data\\kernel\\MEX_HRSC_V09.TI')
      \begintext
