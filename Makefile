VERSION  := $(shell git describe --tags --always --dirty)
BUILD	:= $(shell date +%Y%m%d%H%M%S)

# define image names
REGISTRY := seglh
APP      := basher
DIR := $(shell pwd)

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(VERSION)_rc1.10
IMG_LATEST    := $(IMG):latest

.PHONY: push build tag

push: build tag
	docker push $(IMG_VERSIONED)

build:
	docker buildx build --build-arg IMG_VERSIONED=$(IMG_VERSIONED) --platform linux/amd64 -t $(IMG_VERSIONED) . || \
	docker build  --build-arg IMG_VERSIONED=$(IMG_VERSIONED) -t $(IMG_VERSIONED) .
	docker save $(IMG_VERSIONED) | gzip > $(DIR)/$(REGISTRY)-$(APP):$(VERSION)_$(BUILD).tar.gz

tag:
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)

