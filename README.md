# 粒子フィルターによる自動走行
　
##指定した４つのコーナを通る自動走行　

計算条件  
1 4点のコーナを指定  
2 ロボットは自己位置と走行方向を粒子フィルター推定して走行  
3 目標とするコーナと推定された走行角度より角度を調整  
(https://cloud.githubusercontent.com/assets/20177544/20238236/11a234c0-a92a-11e6-8033-ad6cc5daefc3.png)  
4 コーナの近傍(距離1以内)に達すると次のコーナは向かう  

![Robot trace](https://cloud.githubusercontent.com/assets/20177544/20238207/68f4e714-a929-11e6-81dc-74f64068db33.png)
