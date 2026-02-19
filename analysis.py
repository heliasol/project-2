from abc import ABC, abstractmethod
from parsers import *
from annotations import *
from ontology import *
from hierarchy import *
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt





class GeneAnalyser:
    def __init__(self,annotation_collection: AnnotationCollection,
                 term_collection: TermCollection,
                 hierarchy: OntologyHierarchy):

        self._annotations = annotation_collection
        self._ontology = term_collection
        self._hierarchy = hierarchy


    def _get_ann(self, gene: str) -> list[GeneAnnotation]: #it's private helper
        """
        Internal helper: return all annotations of a gene.
        """
        return self._annotations.get_by_gene_name(gene)
    
    def is_gene_ancestor(self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is not None and a2.term is not None:
                    if self._hierarchy.is_ancestor(a1.term.go_id,a2.term.go_id):
                        return True
        return False

    def is_gene_descendant(self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is not None and a2.term is not None:
                    if self._hierarchy.is_descendant(a1.term.go_id, a2.term.go_id):
                        return True
        return False
    
    def gene_specificity(self, gene: str) -> float | None:
        anns = self._get_ann(gene)
        depths = []
        
        for ann in anns:
            if ann.term is not None:
                ancestors = self._hierarchy._ontology.get_ancestors(ann.term.go_id)
                depths.append(len(ancestors))
        
        if not depths:
            return None

        return sum(depths) / len(depths)


    def genes_functionally_related (self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)

        terms1 = [ann.term.go_id for ann in anns1 if ann.term is not None]
        terms2 = [ann.term.go_id for ann in anns2 if ann.term is not None]


        for t1 in terms1:
            for t2 in terms2:
                if self._hierarchy.is_related(t1, t2):
                    return True
        return False


    def gene_paths(self, gene1: str, gene2: str):
        terms1 = {a.term.go_id for a in self._get_ann(gene1) if a.term}
        terms2 = {a.term.go_id for a in self._get_ann(gene2) if a.term}

        seen = set()
        out = []
        for t1 in terms1:
            for t2 in terms2:
                for p in self._hierarchy.pedigree_paths(t1, t2):
                    if len(p) == 1: #in case both have the same go ids it results in a path length 1
                        continue    #we ignore them
                    tp = tuple(p)
                    if tp not in seen:
                        seen.add(tp) #with seen we avoid returning same paths multiple times
                        out.append(p)
        return out
    
    def shortest_gene_path(self, gene1: str, gene2: str) -> list[str] | None:
        #highly related
        paths = self.gene_paths(gene1, gene2)
        return min(paths, key=len) if paths else None


    def longest_gene_path(self, gene1: str, gene2: str) -> list[str] | None:
       #more complex relationship
       paths = self.gene_paths(gene1, gene2)
       return max(paths, key=len) if paths else None
    
    def MSCA(self, gene1: str, gene2: str) -> str | None:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is None or a2.term is None:
                    continue
                
                msca = self._hierarchy.MSCA(a1.term.go_id, a2.term.go_id)

                if msca:
                    return msca   # return the first common ancestor found

        return None



    

class NumericalAnalysis(ABC):
    def __init__(self, ontology_df : pd.DataFrame, annotation_df : pd.DataFrame):
        self._ontology = ontology_df
        self._annotations = annotation_df

    @property
    @abstractmethod
    def compute(self):
        raise NotImplementedError("You can't instanciate this object -> ABSTRACT!")


class SummaryStatistics(NumericalAnalysis):
    def __init__(self, ontology_df, annotation_df, term: TermCollection):
        super().__init__(ontology_df, annotation_df)
        self.__term = term
        
    @property
    def compute(self):
        onto_df = self._ontology.copy()
        onto_df["n_parents"] = onto_df["parents"].apply(len)
       
        def get_children_info(go_id):
            children = self.__term.get_children(go_id)
            return pd.Series({
                "n_children": len(children),
                "is_leaf": len(children) == 0
            })
        
        
        children_info = onto_df["go_id"].apply(get_children_info)
        onto_df["n_children"] = children_info["n_children"]
        onto_df["is_leaf"] = children_info["is_leaf"]
        leaf_count = onto_df["is_leaf"].sum()
        ann_df = self._annotations.copy()

        experimental_codes = {"EXP","IDA","IPI","IMP","IGI","IEP"}
        ann_df["is_experimental"] = ann_df["evidence"].isin(experimental_codes)
        
        return {
            "namespace counts": onto_df["namespace"].value_counts(),
            "avg parents": onto_df["n_parents"].mean(),
            "avg children": onto_df["n_children"].mean(),
            "leaf_percentage": f"{leaf_count / len(onto_df) * 100:.3f}%",
            "evidence counts": ann_df["evidence"].value_counts(),
            "experimental vs computational": ann_df["is_experimental"].value_counts()
        }
    
    def plots(self):
        summary = self.compute
    
        plt.figure() # Crea una pagina bianca nuova
        summary["namespace counts"].plot(kind="bar", title="Namespace") #ma qua me lo divide per nomi o non me li scrive?
        
        plt.figure() # Crea un'altra pagina bianca nuova
        summary["evidence counts"].plot(kind="bar", title="Evidence")
        
        plt.figure() # Crea la terza pagina bianca
        summary["experimental vs computational"].plot(kind="bar", title="Exp vs Comp")
        
        plt.show() # Mostra tutto


'''principle: the more go_term 2 gene share the more they are related'''
class GeneSimilarityAnalysis(NumericalAnalysis):
    def __init__(self, ontology_df, annotation_df):
        super().__init__(ontology_df, annotation_df)
        self.__sim = None  #so it doesn't get built from scratch every time
        self.__gene2terms =None
    @property
    def compute(self):
        """
        Returns gene-gene Jaccard similarity matrix (limited to top 500 most annotated genes)
        This makes computation fast (~5 seconds) instead of hours
        """
        if self.__sim is not None:
            return self.__sim     #if it is already built it just returns it

        # gene × term table (binary)   1 se gena ha quel GO c'è - 0 se non lo ha
        table = pd.crosstab(
            self._annotations["gene_id"],
            self._annotations["go_id"]
        )
        
        # most 500 annotated genes
        # counts annotation per gene -> table with decreasing numbers top to bottom
        gene_counts = table.sum(axis=1)
        # most 500 annotated genes
        top_genes = gene_counts.nlargest(500).index
        # filters the table just to keep that 500 genes
        table = table.loc[top_genes]
        
        # convert to numpy, because it's easier to work with matricial and logical operations
        M = table.values.astype(bool)  # boolean faster

        n = M.shape[0]  # n = 500 
        sim = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                intersection = np.logical_and(M[i], M[j]).sum()
                union = np.logical_or(M[i], M[j]).sum()

                sim[i, j] = intersection / union if union else 0  #jaccard similarity

        self.__sim = pd.DataFrame(sim, index=table.index, columns=table.index)
        return self.__sim

    def compare2genes(self, gene1, gene2):
        if self.__gene2terms is None:
            self.__gene2terms = (
                self._annotations.groupby('gene_id')['go_id'].apply(set).to_dict()
            )

        gene_1 = self.__gene2terms.get(gene1)
        gene_2 = self.__gene2terms.get(gene2)

        if not gene_1 or not gene_2:
            return 0.0

        inter = len(gene_1 & gene_2)
        union = len(gene_1 | gene_2)
        return round(inter / union, 3) if union else 0.0
        



    






