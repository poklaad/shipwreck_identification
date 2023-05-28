# shipwreck_identification
Программная реализация распозавания типовых аварийных состояний судна

## Описание
Программная реализация представлена в разделах **code** и **neural_network**.

### code 
Содержит в себе программную реализацию модели поведения корабля в результате развития аварии. 
**Dynamic_recovery.m** - функция восстанавливающего момента судна.
**Static_oscillation_with_phase_portrait.m** - график угла крена судна и фазового портрета при статичном аварийном состоянии.
**Oscillation_with_phase_portrait.m** - график угла крена судна и фазового портрета при развитии аварии.
**data_modeling.m** - реализация модели для построения набора данных для обучения нейронных сетей.

### neural_network 
Содержит в себе файлы из Google Colab, в которых обучались нейронные сети.
**multilayer_perceptron.ipynb** - нейронная сеть на основе многослойного перцептрона.
**Kohonen_network.ipynb** - нейронная сеть на основе сети Кохонена.
**hybrid_network.ipynb** - нейронная сеть на основе гибридной сети.

### dataset 
Содержит в себе данные для обучения нейронных сетей.
**balanced_autocor.txt** - данные о функции автокорреляции. На них обучаются нейронные сети.
**balanced_wreck_type.txt** - метки, используемые при обучении. Представляют собой типы затопления, соответствующие функции автокорреляции

### waving
Содержит в себе данные о различных типах волнения, типичных для одного из районов Баренцева моря.
**ANG4.DAT** - волновой склон 4 балла
**ANG5.DAT** - волновой склон 5 баллов
**ANG6.DAT** - волновой склон 6 баллов
**ANG7.DAT** - волновой склон 7 баллов
**ANG8.DAT** - волновой склон 8 баллов
**anglM.DAT** - ветровое волнение
**anglS.DAT** - зыбь
**anglWW.DAT** - смешанное волнение
