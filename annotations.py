from ontology import *

class GeneAnnotation:
        branch_map = {"P": "Biological Process", "F": "Molecular Function", "C": "Cellular Component"}  

        def __init__(self, gene_id: str, 
                     gene_name: str,     
                     go_id: str, 
                     qualifier: str | None = None,  
                     aspect: str | None = None, 
                     evidence: str | None = None, 
                     molecule: str = "protein" ):    
            
            self.__gene_id = gene_id
            self.__gene_name = gene_name
            self.__go_id = go_id
            self.__qualifier = qualifier
            self.__aspect = aspect
            self.__evidence = evidence
            self.__molecule = molecule
            self.__branch = self.branch_map.get(aspect, "Unknown")

            self.__term: Term | None = None     

        @property
        def gene_name(self) -> str:
            return self.__gene_name

        @property
        def go_id(self) -> str:
            return self.__go_id

        @property
        def qualifier(self) -> str | None:
            return self.__qualifier

        @property
        def aspect(self) -> str | None:
            return self.__aspect

        @property
        def branch(self) -> str:
            return self.__branch

        @property
        def evidence(self) -> str | None:
            return self.__evidence

        @property
        def molecule(self) -> str:
            return self.__molecule

        @property 
        def term(self) -> Term | None:
            return self.__term

        def link_term(self, term_collection: "TermCollection") -> None:
            term = term_collection.get_term(self.go_id)
            if term != None:
                self.__term = term


        def __repr__(self) -> str:
            return f' gene name: {self.gene_name} \n go id: {self.go_id} \n branch: {self.branch} \n evidence: {self.evidence} \n molecule: {self.molecule} \n'


class AnnotationCollection: 

    def __init__(self) -> None:
        # internal mutable storage
        self._annotations: list[GeneAnnotation] = []

    @property
    def annotations(self) -> list[GeneAnnotation]:
        return self._annotations.copy()

    def __len__(self) -> int:
        return len(self._annotations)

    def __iter__(self):
        return iter(self._annotations)

    def __repr__(self) -> str:
        return f"<AnnotationCollection: {len(self._annotations)} annotations>"

    def add_annotation(self, annotation: GeneAnnotation) -> None:
        self._annotations.append(annotation)

    def link_terms(self, term_collection: "TermCollection") -> None:
        for ann in self.annotations:
            ann.link_term(term_collection)

    def get_by_gene_id(self, gene_id: str) -> list[GeneAnnotation]: 
        return [ann for ann in self._annotations if ann.gene_id == gene_id]

    def get_by_gene_name(self, gene_name: str) -> list[GeneAnnotation]:  
        return [ann for ann in self._annotations if ann.gene_name == gene_name]

    def get_by_term(self, go_id: str) -> list[GeneAnnotation]: 
        return [ann for ann in self._annotations if ann.go_id == go_id]

    def get_by_aspect(self, aspect: str) -> list[GeneAnnotation]: 
        return [ann for ann in self._annotations if ann.aspect == aspect]

    def get_by_evidence(self, evidence: str) -> list[GeneAnnotation]: 
        return [ann for ann in self._annotations if ann.evidence == evidence]

    def __repr__(self):
        return f"<AnnotationCollection: {len(self.annotations)} annotations>"

