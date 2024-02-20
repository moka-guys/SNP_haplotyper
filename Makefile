VERSION  := $(shell git describe --tags --always --dirty)

# define image names
REGISTRY := seglh
APP      := basher

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(VERSION)
IMG_LATEST    := $(IMG):latest

.PHONY: push build tag

push: build tag
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

build:
	docker buildx build --build-arg IMG_VERSIONED=$(IMG_VERSIONED) --platform linux/amd64 -t $(IMG_VERSIONED) . || \
	docker build  --build-arg IMG_VERSIONED=$(IMG_VERSIONED) -t $(IMG_VERSIONED) .

tag:
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)

