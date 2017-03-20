const static double kinv6[625]={0.47606736897589650250108005230538473,-0.201999054135754018879300868242938103,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.04120416580489373123761450639677008,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,0.047489769608872611307421749937661943,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,0.0259329067334364259339882294506009867,-0.0212802077523094021808359624598092963,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.208560981968879116524638471810649902,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,-0.201999054135754018879300868242938103,0.50938573318431102969008756056795776,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.043888547718481721701661965177656815,-0.171617696998178203155779765713693344,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.0232213372836719097279322344720453274,0.04512844678562585311935256370572501,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,-0.0226798740595273063381248993584034593,0.0156197644369657136399942587003296616,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.30174861502680050916594581239208528,-0.029076871581563580834947425304830934,-0.055921078150567034253562327924698337,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.46942845198569621631777061736028368,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.0086024374158155004863305786489671871,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,-0.209392499754775496425226410819288084,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,-0.0226798740595273063381248993584034593,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.077669658481918559119834869425197733,-0.126186637712880924912638402676012118,-0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,-0.208705124825054220766464055078349331,0.47606736897589650250108005230538473,0.0156197644369657136399942587003296616,0.04120416580489373123761450639677008,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.047489769608872611307421749937661943,-0.208705124825054220766464055078349331,-0.0212802077523094021808359624598092963,0.0156197644369657136399942587003296616,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.034427444653135554543702713760963884,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,-0.201999054135754018879300868242938103,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,0.0156197644369657136399942587003296616,0.50938573318431102969008756056795776,-0.171617696998178203155779765713693344,-0.201999054135754018879300868242938103,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,0.0156197644369657136399942587003296616,0.0232213372836719097279322344720453274,-0.0226798740595273063381248993584034593,-0.201999054135754018879300868242938103,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,0.0156197644369657136399942587003296616,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,0.30174861502680050916594581239208528,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,-0.055921078150567034253562327924698337,0.04120416580489373123761450639677008,-0.171617696998178203155779765713693344,0.0086024374158155004863305786489671871,0.04120416580489373123761450639677008,-0.171617696998178203155779765713693344,0.65082689318597309796008104638122202,0.04120416580489373123761450639677008,-0.171617696998178203155779765713693344,0.0086024374158155004863305786489671871,0.04120416580489373123761450639677008,-0.031439860068336528554438658839402338,0.0086024374158155004863305786489671871,0.04120416580489373123761450639677008,-0.171617696998178203155779765713693344,0.0086024374158155004863305786489671871,0.04120416580489373123761450639677008,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,0.97575060417554501098578771906700104,0.128773269508504881924166907389650801,-0.128773269508504881924166907389650801,0.128773269508504881924166907389650801,-0.128773269508504881924166907389650801,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.04120416580489373123761450639677008,0.47606736897589650250108005230538473,-0.201999054135754018879300868242938103,0.0259329067334364259339882294506009867,-0.0212802077523094021808359624598092963,0.047489769608872611307421749937661943,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,0.034427444653135554543702713760963884,-0.094867545488743565676208572421196756,0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,-0.079265991826999996304727185628489262,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.043888547718481721701661965177656815,-0.171617696998178203155779765713693344,-0.201999054135754018879300868242938103,0.50938573318431102969008756056795776,-0.0226798740595273063381248993584034593,0.0156197644369657136399942587003296616,0.0232213372836719097279322344720453274,0.04512844678562585311935256370572501,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,0.21675066529466989407743605916255601,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,0.055921078150567034253562327924698337,0.029076871581563580834947425304830934,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,0.04512844678562585311935256370572501,0.0086024374158155004863305786489671871,0.0259329067334364259339882294506009867,-0.0226798740595273063381248993584034593,0.46942845198569621631777061736028368,-0.208705124825054220766464055078349331,-0.209392499754775496425226410819288084,0.048348773637094730318427987761544205,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.077669658481918559119834869425197733,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,-0.126186637712880924912638402676012118,-0.090815707806348549388871629006537109,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,0.04120416580489373123761450639677008,-0.0212802077523094021808359624598092963,0.0156197644369657136399942587003296616,-0.208705124825054220766464055078349331,0.47606736897589650250108005230538473,0.047489769608872611307421749937661943,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.034427444653135554543702713760963884,0.079265991826999996304727185628489262,-0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.094867545488743565676208572421196756,0.047489769608872611307421749937661943,0.0232213372836719097279322344720453274,-0.209392499754775496425226410819288084,0.047489769608872611307421749937661943,0.0232213372836719097279322344720453274,-0.031439860068336528554438658839402338,0.047489769608872611307421749937661943,0.0232213372836719097279322344720453274,-0.209392499754775496425226410819288084,0.047489769608872611307421749937661943,0.46822666098313243800020237145874682,-0.209392499754775496425226410819288084,0.047489769608872611307421749937661943,0.0232213372836719097279322344720453274,-0.209392499754775496425226410819288084,0.047489769608872611307421749937661943,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,-0.035227500154945551450248342065137204,0.097843568762155153181163906469981746,-0.097843568762155153181163906469981746,0.097843568762155153181163906469981746,-0.097843568762155153181163906469981746,0.0259329067334364259339882294506009867,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,-0.208705124825054220766464055078349331,-0.0226798740595273063381248993584034593,0.0086024374158155004863305786489671871,0.0259329067334364259339882294506009867,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,-0.208705124825054220766464055078349331,-0.209392499754775496425226410819288084,0.46942845198569621631777061736028368,0.0259329067334364259339882294506009867,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,-0.208705124825054220766464055078349331,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,-0.139332687037310915181675162257351494,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,0.126186637712880924912638402676012118,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,0.0259329067334364259339882294506009867,-0.0212802077523094021808359624598092963,-0.201999054135754018879300868242938103,0.04120416580489373123761450639677008,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,0.047489769608872611307421749937661943,0.0259329067334364259339882294506009867,0.47606736897589650250108005230538473,-0.201999054135754018879300868242938103,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.034427444653135554543702713760963884,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,-0.0226798740595273063381248993584034593,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,-0.171617696998178203155779765713693344,0.0156197644369657136399942587003296616,0.043888547718481721701661965177656815,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.0232213372836719097279322344720453274,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.50938573318431102969008756056795776,0.04512844678562585311935256370572501,-0.201999054135754018879300868242938103,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.21675066529466989407743605916255601,0.055921078150567034253562327924698337,0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,0.0259329067334364259339882294506009867,-0.0226798740595273063381248993584034593,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,0.04512844678562585311935256370572501,0.0086024374158155004863305786489671871,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.048348773637094730318427987761544205,0.0259329067334364259339882294506009867,-0.209392499754775496425226410819288084,0.048348773637094730318427987761544205,-0.208705124825054220766464055078349331,0.04512844678562585311935256370572501,0.46942845198569621631777061736028368,-0.208705124825054220766464055078349331,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.139332687037310915181675162257351494,0.090815707806348549388871629006537109,0.126186637712880924912638402676012118,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,-0.0212802077523094021808359624598092963,0.0156197644369657136399942587003296616,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.0156197644369657136399942587003296616,0.04120416580489373123761450639677008,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,0.0259329067334364259339882294506009867,0.047056700500673124240131459797551201,0.047489769608872611307421749937661943,-0.208705124825054220766464055078349331,0.047056700500673124240131459797551201,-0.201999054135754018879300868242938103,-0.208705124825054220766464055078349331,0.47606736897589650250108005230538473,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.139706092662608007437233044288722134,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,-0.0092273327027093340061026585412559696,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,-0.053637747404387036538383169415883882,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,0.0057398078290956923747206232859137222,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,0.061995248417598442138174881428363804,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.090323753809315826966155354758880816,-0.3359902541872321754844773217002408,-0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,0.0241941425336893278105239728891133004,0.061995248417598442138174881428363804,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,-0.036072620747944560269137923139530078,-0.053637747404387036538383169415883882,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,0.0057398078290956923747206232859137222,0.0103846033740640180498128650637113438,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,-0.0092273327027093340061026585412559696,-0.0208965840633061976516446304038626469,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.44005638392646861361426405374743028,0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,-0.036072620747944560269137923139530078,-0.053637747404387036538383169415883882,0.0241941425336893278105239728891133004,0.061995248417598442138174881428363804,-0.0092273327027093340061026585412559696,-0.0208965840633061976516446304038626469,0.0057398078290956923747206232859137222,0.0103846033740640180498128650637113438,-0.0208965840633061976516446304038626469,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,0.0241941425336893278105239728891133004,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.090323753809315826966155354758880816,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,-0.0137423759299206111636313772883086637,0.0241941425336893278105239728891133004,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,-0.0208965840633061976516446304038626469,0.061995248417598442138174881428363804,-0.053637747404387036538383169415883882,0.0241941425336893278105239728891133004,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,-0.0208965840633061976516446304038626469,0.0057398078290956923747206232859137222,-0.0092273327027093340061026585412559696,0.0241941425336893278105239728891133004,-0.036072620747944560269137923139530078,0.0103846033740640180498128650637113438,-0.0208965840633061976516446304038626469,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.44005638392646861361426405374743028,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,0.208560981968879116524638471810649902,0.30174861502680050916594581239208528,0.077669658481918559119834869425197733,0.034427444653135554543702713760963884,0.30174861502680050916594581239208528,0.97575060417554501098578771906700104,0.034427444653135554543702713760963884,0.21675066529466989407743605916255601,0.077669658481918559119834869425197733,0.034427444653135554543702713760963884,-0.035227500154945551450248342065137204,-0.139332687037310915181675162257351494,0.034427444653135554543702713760963884,0.21675066529466989407743605916255601,-0.139332687037310915181675162257351494,-0.139706092662608007437233044288722134,-0.090323753809315826966155354758880816,-0.44005638392646861361426405374743028,-0.090323753809315826966155354758880816,-0.44005638392646861361426405374743028,6.4474518494168764557121225501868598,-0.55664515567372376038704731943137047,-3.0660029894550230953962415843578202,-0.55664515567372376038704731943137047,-3.0660029894550230953962415843578202,-0.094867545488743565676208572421196756,-0.029076871581563580834947425304830934,-0.126186637712880924912638402676012118,-0.094867545488743565676208572421196756,-0.029076871581563580834947425304830934,0.128773269508504881924166907389650801,-0.094867545488743565676208572421196756,-0.029076871581563580834947425304830934,0.090815707806348549388871629006537109,0.079265991826999996304727185628489262,0.097843568762155153181163906469981746,0.090815707806348549388871629006537109,0.079265991826999996304727185628489262,0.055921078150567034253562327924698337,0.090815707806348549388871629006537109,0.079265991826999996304727185628489262,-0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.55664515567372376038704731943137047,3.824188374201590623804406888784741,-0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,-0.079265991826999996304727185628489262,-0.055921078150567034253562327924698337,-0.090815707806348549388871629006537109,-0.079265991826999996304727185628489262,0.029076871581563580834947425304830934,-0.128773269508504881924166907389650801,0.094867545488743565676208572421196756,0.029076871581563580834947425304830934,-0.090815707806348549388871629006537109,-0.079265991826999996304727185628489262,-0.097843568762155153181163906469981746,-0.090815707806348549388871629006537109,0.094867545488743565676208572421196756,0.029076871581563580834947425304830934,0.126186637712880924912638402676012118,0.094867545488743565676208572421196756,-0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,-0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-3.0660029894550230953962415843578202,-0.201540229072843768021117984995550364,3.824188374201590623804406888784741,-0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,-0.094867545488743565676208572421196756,-0.029076871581563580834947425304830934,0.090815707806348549388871629006537109,0.079265991826999996304727185628489262,-0.029076871581563580834947425304830934,0.128773269508504881924166907389650801,0.079265991826999996304727185628489262,0.055921078150567034253562327924698337,-0.126186637712880924912638402676012118,-0.094867545488743565676208572421196756,0.097843568762155153181163906469981746,0.090815707806348549388871629006537109,-0.094867545488743565676208572421196756,-0.029076871581563580834947425304830934,0.090815707806348549388871629006537109,0.079265991826999996304727185628489262,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,-0.55664515567372376038704731943137047,0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,3.824188374201590623804406888784741,-0.201540229072843768021117984995550364,-0.079265991826999996304727185628489262,0.029076871581563580834947425304830934,-0.090815707806348549388871629006537109,0.094867545488743565676208572421196756,-0.055921078150567034253562327924698337,-0.128773269508504881924166907389650801,-0.079265991826999996304727185628489262,0.029076871581563580834947425304830934,-0.090815707806348549388871629006537109,0.094867545488743565676208572421196756,-0.097843568762155153181163906469981746,0.126186637712880924912638402676012118,-0.079265991826999996304727185628489262,0.029076871581563580834947425304830934,-0.090815707806348549388871629006537109,0.094867545488743565676208572421196756,-0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,-3.0660029894550230953962415843578202,-0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,3.824188374201590623804406888784741};
