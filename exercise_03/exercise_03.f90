PROGRAM evolution

	IMPLICIT NONE	!// Prevent the use of implicit declaration
	
	INTEGER							:: nSpecies = 64
	REAL, DIMENSION(64) 			:: chain
	INTEGER							:: tPeriod = 40
	INTEGER, DIMENSION(20)			:: mAct				!// Where the mutation activity happens
	INTEGER							:: t = 1
	INTEGER							:: i
	REAL							:: minBarrier
	
	DO i = 1, nSpecies
		CALL RANDOM_NUMBER(chain(i))
	END DO
	
	WRITE(*,*) SUM(chain)/(MAX(1,SIZE(chain)))
	
	DO WHILE (t <= tPeriod)
		minBarrier = MINVAL(chain)
			DO i = 1, nSpecies
				IF (chain(i) == minBarrier) THEN
					mAct(t) = i
					CALL RANDOM_NUMBER(chain(i))
		
					IF (i == 1) THEN
						CALL RANDOM_NUMBER(chain(nSpecies))
					ELSE
						CALL RANDOM_NUMBER(chain(i-1))
					END IF
					
					IF (i == nSpecies) THEN
						CALL RANDOM_NUMBER(chain(1))
					ELSE
						CALL RANDOM_NUMBER(chain(i+1))
					END IF
					
					EXIT
				END IF
			END DO
		t = t + 1
	END DO
	
	WRITE(*,*) SUM(chain)/(MAX(1,SIZE(chain)))
	
END PROGRAM evolution
