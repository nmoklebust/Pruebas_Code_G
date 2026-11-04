"""Audita la selección descriptiva de los modelos mostrados en las mezclas nuevas.

Uso, desde la raíz del repositorio:
    python seleccion_modelos_mezclas/calcular_metricas.py
Requiere numpy y pandas; lee figuras_poster.zip sin extraerlo ni modificarlo.
"""
from pathlib import Path
import argparse, ast, hashlib, json, zipfile
import numpy as np
import pandas as pd

NOMBRES={'Linear_Regression':'OLS','PLS_Regression':'PLS','Ridge_Regression':'Ridge',
         'SVR_RBF':'SVR (RBF)','PCA___GPR':'GPR-PCA'}
OBJETIVOS={'NaCl':'NaCl nominal','KCl':'KCl nominal','PBS':'Cloruro total nominal (NaCl + KCl)'}
ORDEN={'NaCl':0,'PBS':1,'KCl':2}

def calcular(paquete,salida):
    paquete,salida=Path(paquete),Path(salida);salida.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(paquete) as z:
        prefijo='figuras_poster/'
        origen=prefijo+'salidas/tablas/predicciones_mezclas.csv'
        sha_original=hashlib.sha256(z.read(origen)).hexdigest()
        with z.open(origen) as f:p=pd.read_csv(f)
        modelos=json.loads(z.read(prefijo+'salidas/tablas/modelos_inferencia.json'))
        refs=json.loads(z.read(prefijo+'salidas/tablas/referencias_porcentaje.json'))
        madres={r['solucion']:r['madre_mM'] for r in refs['referencias']}
        fuente=z.read(prefijo+'figuras_poster.py').decode('utf-8')
        tree=ast.parse(fuente)
        node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SELECCIONADOS' for t in n.targets))
        mostrados=ast.literal_eval(node.value)
        clave=pd.read_csv(z.open(prefijo+'salidas/tablas/clave_mezclas.csv'))
    if len(p)!=2700 or p.duplicated(['familia','modelo','row_index']).any():
        raise ValueError('La tabla no contiene las 2700 predicciones individuales únicas esperadas')
    modos={(m['solucion'],m['modelo']):m['modo'] for m in modelos}
    filas=[];por_mezcla=[]
    for familia,grupo_familia in p.groupby('familia'):
        if familia not in OBJETIVOS:raise ValueError('Objetivo no reconocido')
        referencia=None
        for modelo,g in grupo_familia.groupby('modelo'):
            g=g.sort_values('row_index')
            identidad=g[['row_index','sample_id','mezcla','valor_real_mM']].reset_index(drop=True)
            if referencia is None:referencia=identidad
            else:pd.testing.assert_frame_equal(identidad,referencia)
            if len(g)!=300 or not g.groupby('mezcla').size().eq(60).all():
                raise ValueError('Las mezclas o recuentos no coinciden con la evaluación de referencia')
            esperado=(g.nominal_nacl_mM if familia=='NaCl' else g.nominal_kcl_mM if familia=='KCl'
                      else g.nominal_nacl_mM+g.nominal_kcl_mM)
            np.testing.assert_allclose(g.valor_real_mM,esperado,atol=1e-10,rtol=0)
            real=g.valor_real_mM.to_numpy();pred=g.valor_predicho_mM.to_numpy();error=pred-real
            np.testing.assert_allclose(error,g.error_mM,atol=1e-10,rtol=0)
            mae=float(np.mean(np.abs(error)));rmse=float(np.sqrt(np.mean(error**2)))
            sst=float(np.sum((real-real.mean())**2))
            filas.append(dict(familia=familia,objetivo=OBJETIVOS[familia],modelo=NOMBRES[modelo],
                identificador_modelo=modelo,modo='Optimizado' if modos[familia,modelo]=='search' else 'Base',
                n=300,n_mezclas=5,n_por_mezcla=60,MAE_mM=mae,RMSE_mM=rmse,
                sesgo_mM=float(error.mean()),R2=float(1-np.sum(error**2)/sst),
                madre_mM=madres[familia],MAE_pp=100*mae/madres[familia],
                mostrado_en_nominal_vs_predicha=mostrados[familia]==modelo))
            for mezcla,h in g.groupby('mezcla'):
                e=h.valor_predicho_mM-h.valor_real_mM
                por_mezcla.append(dict(familia=familia,modelo=NOMBRES[modelo],identificador_modelo=modelo,
                    mezcla=mezcla,n=len(h),nominal_mM=float(h.valor_real_mM.iloc[0]),
                    MAE_mM=float(e.abs().mean()),sesgo_mM=float(e.mean())))
    globales=pd.DataFrame(filas)
    globales['rango_MAE_en_objetivo']=globales.groupby('familia').MAE_mM.rank(method='min',ascending=True).astype(int)
    globales['_orden']=globales.familia.map(ORDEN)
    globales=globales.sort_values(['_orden','rango_MAE_en_objetivo']).drop(columns='_orden')
    if not (globales.mostrado_en_nominal_vs_predicha==(globales.rango_MAE_en_objetivo==1)).all():
        raise ValueError('La selección fijada en la figura no coincide con el mínimo MAE calculado')
    detalle=pd.DataFrame(por_mezcla)
    ganadores=detalle.loc[detalle.groupby(['familia','mezcla']).MAE_mM.idxmin()].copy()
    seleccion=[]
    for familia,g in globales.groupby('familia',sort=False):
        primero,segundo=g.iloc[0],g.iloc[1]
        macro=detalle[(detalle.familia==familia)&(detalle.modelo==primero.modelo)].MAE_mM.mean()
        np.testing.assert_allclose(macro,primero.MAE_mM,atol=1e-10,rtol=0)
        seleccion.append(dict(familia=familia,objetivo=primero.objetivo,modelo_mostrado=primero.modelo,
            modo=primero.modo,MAE_global_mM=primero.MAE_mM,RMSE_global_mM=primero.RMSE_mM,
            segundo_modelo=segundo.modelo,MAE_segundo_mM=segundo.MAE_mM,
            diferencia_MAE_mM=segundo.MAE_mM-primero.MAE_mM,
            reduccion_MAE_vs_segundo_pct=100*(segundo.MAE_mM-primero.MAE_mM)/segundo.MAE_mM,
            mezclas_con_menor_MAE=int(((ganadores.familia==familia)&(ganadores.modelo==primero.modelo)).sum()),
            n=300,seleccion='descriptiva posterior a la evaluación de estas mezclas'))
    for nombre,df in [('metricas_globales_nueve_modelos.csv',globales),
                      ('MAE_por_mezcla_nueve_modelos.csv',detalle),
                      ('ganador_por_mezcla_y_objetivo.csv',ganadores),
                      ('seleccion_por_objetivo.csv',pd.DataFrame(seleccion)),
                      ('concentraciones_nominales_mezclas.csv',clave)]:
        df.to_csv(salida/nombre,index=False,float_format='%.17g')
    # Resumen compacto solicitado para discutir MAE, RMSE y R² de las mezclas nuevas.
    columnas=['familia','objetivo','modelo','modo','n','MAE_mM','RMSE_mM','R2','sesgo_mM',
              'mostrado_en_nominal_vs_predicha']
    globales[columnas].rename(columns={'familia':'solucion_entrenamiento','modo':'configuracion',
        'n':'n_observaciones','mostrado_en_nominal_vs_predicha':'modelo_mostrado_en_figura'}).to_csv(
        salida/'metricas_MAE_RMSE_R2_mezclas_nuevas.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
    info=dict(paquete=paquete.name,tabla_entrada=origen,
              sha256_predicciones=sha_original,
              n_predicciones=2700,n_modelos=9,n_por_modelo=300,mezclas_por_modelo=5,
              observaciones_por_mezcla_y_modelo=60,objetivos_verificados=OBJETIVOS,
              criterio='MAE global mínimo dentro de cada objetivo; selección descriptiva para visualización',
              formula='MAE = mean(abs(predicho_mM - nominal_mM)) sobre 300 observaciones',
              equivalencia='Los cinco MAE por mezcla tienen el mismo peso porque cada mezcla aporta 60 observaciones',
              lista_fijada_en_codigo=mostrados,lista_coincide_con_minimos_MAE=True,
              unidades='mM; MAE_pp = 100*MAE_mM/madre_mM',
              limitacion='Elección posterior a observar el desempeño de estas mezclas, no una selección independiente previa; observaciones repetidas descriptivas')
    (salida/'procedencia_y_calculo.json').write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n')
    print(globales[['familia','modelo','modo','MAE_mM','RMSE_mM','sesgo_mM','R2','mostrado_en_nominal_vs_predicha']].to_string(index=False))
    print(pd.DataFrame(seleccion).to_string(index=False))
    return globales

if __name__=='__main__':
    raiz=Path(__file__).resolve().parent.parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paquete',type=Path,default=raiz/'figuras_poster.zip')
    parser.add_argument('--salida',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();calcular(args.paquete,args.salida)
