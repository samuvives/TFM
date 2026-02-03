
    # ===================================================
    # 3.10 RANKING GLOBAL
    # ===================================================
    print(f"[3.10] Calculando ranking global para K={k}")
    all_feats = []
    for view, W in weights.items():
        imp = (W**2).sum(axis=0)
        df_imp = imp.to_frame(name='GlobalScore')
        df_imp['View'] = view
        df_imp['PctContrib'] = (df_imp['GlobalScore'] / df_imp['GlobalScore'].sum()) * 100
        all_feats.append(df_imp)

    # Definimos df_global uniendo todos los datos
    df_global = (pd
        .concat(all_feats)
        .sort_values(by='GlobalScore', ascending=False))
    print("todo okey en 3.10")

    # ===================================================
    # 3.11 LIDERES POR VISTA (TOP 3)
    # ===================================================
    print(f"[3.11] Extrayendo los Top 3 líderes de cada vista")

    # Reseteamos el índice para que el nombre del gen/metabolito sea una columna llamada 'Feature'
    df_global_reset = df_global.reset_index().rename(columns={'index': 'Feature'})

    # Extraemos el top 3 de cada grupo 'View'
    top_leaders = (
        df_global_reset
        .groupby('View')
        .apply(lambda x: x.sort_values('GlobalScore', ascending=False).head(3))
        .reset_index(drop=True)
    )

    top_leaders.to_csv(os.path.join(outdir, "top_3_leaders_per_view.csv"), index=False)
    print("todo okey en 3.11")

    # ===================================================
    # 3.12 HEATMAP DE LÍDERES (INTEGRADO)
    # ===================================================
    print("[3.12] Generando Heatmap de líderes")
    list_w = []
    for _, row in top_leaders.iterrows():
        v = row['View']
        f = row['Feature']

        # Accedemos a los pesos originales usando la vista y el nombre de la variable
        ws = weights[v][f]
        ws.name = f"{v}: {f}"
        list_w.append(ws)

    plt.figure(figsize=(12, 10))
    sns.heatmap(pd.concat(list_w, axis=1).T, cmap="RdBu_r", center=0)
    plt.title(f"Integrated Leaders Heatmap (K={k})")
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "integrated_leaders_heatmap.png"))
    plt.close()
    print("todo okey en 3.12")

    # ===================================================
    # 3.13 LOADING PLOTS (ENFOQUE DE LÍDERES) - BLINDADO
    # ===================================================
    print(f"[3.13] Generando Loading Plots para los líderes de cada vista (K={k})")

    # Calculamos el número de vistas para el layout
    n_views = len(weights)
    fig, axes = plt.subplots(n_views, 1, figsize=(10, 4 * n_views), sharex=False)

    # Ajuste por si solo hay una vista (axes no sería una lista)
    if n_views == 1:
        axes = [axes]

    for i, (view, W) in enumerate(weights.items()):
        # Buscamos el líder #1 de esta vista en nuestro DataFrame top_leaders
        try:
            # Filtramos el top_leaders para esta vista y tomamos el primero
            leader_row = top_leaders[top_leaders['View'] == view].iloc[0]
            leader_name = leader_row['Feature']

            # Extraemos los pesos del líder (Series: Factores -> Pesos)
            leader_loadings = W[leader_name]

            # Graficamos
            sns.barplot(x=leader_loadings.index, y=leader_loadings.values, ax=axes[i], palette="vlag")

            axes[i].set_title(f"Líder de {view}: {leader_name}", fontsize=12, fontweight='bold')
            axes[i].set_ylabel("Weight (W)")
            axes[i].axhline(0, color='black', linewidth=0.8)
            axes[i].tick_params(axis='x', rotation=45)
            print("todo okey en 3.13")

        except Exception as e:
            print(f"   [!] Error procesando líder para la vista {view}: {e}")
            axes[i].set_title(f"Error en vista {view}")

    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "leaders_loading_focus.png"), dpi=300)
    plt.close()

